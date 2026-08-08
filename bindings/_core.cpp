#include "hammock/abstract_sketch.hpp"
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/processing_modes.hpp"
#include "hammock/xxhash.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <exception>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace py = pybind11;

namespace {

SubBMethod parse_subB_method(const std::string& s) {
    if (s == "hash-threshold") return SubBMethod::HashThreshold;
    if (s == "mixed-stride")   return SubBMethod::MixedStride;
    if (s == "single-hash") return SubBMethod::SingleHash;
    throw std::invalid_argument(
        "subB_method must be one of 'hash-threshold', 'mixed-stride', "
        "'single-hash' (got '" + s + "')");
}

// Build an HLLSketch from one BED file in the chosen mode.
HLLSketch sketch_bed_file_hll(const std::string& path,
                              const std::string& mode,
                              size_t precision,
                              const std::string& separator,
                              double sub_a,
                              double sub_b,
                              double exp_a,
                              const std::string& subB_method,
                              uint64_t seed,
                              uint32_t gate_seed,
                              bool verbose) {
    const SubBMethod method = parse_subB_method(subB_method);
    HLLSketch sketch(precision);
    if (mode == "A") {
        process_bed_file_mode_a(path, sketch, seed, separator, verbose);
    } else if (mode == "B") {
        process_bed_file_mode_b(path, sketch, seed, separator, sub_b, method,
                                gate_seed, verbose);
    } else if (mode == "C") {
        process_bed_file_mode_c(path, sketch, seed, separator, sub_a, sub_b, exp_a,
                                method, gate_seed, verbose);
    } else {
        throw std::invalid_argument("Mode must be one of A, B, C (got '" + mode + "')");
    }
    return sketch;
}


// An exception may not propagate out of an OpenMP structured block: escaping a
// `#pragma omp parallel for` calls std::terminate rather than unwinding, so a
// mismatched-precision throw inside the pair loops used to SIGABRT the whole
// interpreter instead of surfacing as a Python RuntimeError. Validating the
// operands up front keeps the common failure out of the parallel region
// entirely; OmpError below catches whatever else might still throw.
void require_uniform_precision(const std::vector<const HLLSketch*>& a,
                               const std::vector<const HLLSketch*>& b) {
    // Null scan first, and *before* the empty short-circuit below. The lists
    // hold borrowed pointers now, and pybind's pointer caster turns `None` into
    // a null rather than refusing it. `.noconvert()` on the bindings already
    // makes None fail to load, so reaching here with a null should be
    // impossible -- but a null that got through would be dereferenced with the
    // GIL released inside an OpenMP region, i.e. a SIGSEGV with no traceback
    // and no chance for OmpError to latch it. Note the check cannot live after
    // the short-circuit: `(a=[s, None], b=[])` returns early there while the
    // cardinality prepass still walks every element of `a`.
    for (const auto& v : {std::cref(a), std::cref(b)}) {
        for (const HLLSketch* s : v.get()) {
            if (!s) throw std::invalid_argument("sketch list contains a null entry");
        }
    }
    if (a.empty() || b.empty()) return;
    const size_t p = a[0]->precision();
    const size_t h = a[0]->hash_size_bits();
    for (const auto& v : {std::cref(a), std::cref(b)}) {
        for (const HLLSketch* s : v.get()) {
            if (s->precision() != p || s->hash_size_bits() != h) {
                throw std::invalid_argument(
                    "all sketches must share precision and hash size; got p=" +
                    std::to_string(s->precision()) + " alongside p=" +
                    std::to_string(p));
            }
        }
    }
}

// Team size for the pairwise loops. `threads <= 0` means "whatever OpenMP would
// have picked", so passing omp_get_max_threads() back in is a no-op.
//
// Without this the pairwise phase ignored --threads entirely: omp_set_num_threads
// appears only in hammock_cli.cpp, so `hammock --threads 4` inside a 4-CPU cgroup
// still ran this loop on every core on the node. Applied as a num_threads()
// clause rather than omp_set_num_threads() so it cannot leak into unrelated
// regions -- notably the sketching loops in processing_modes.cpp, which are
// nested inside a Python thread pool and would compound rather than clamp.
int pairwise_team_size(int threads) {
#ifdef _OPENMP
    return (threads > 0) ? threads : omp_get_max_threads();
#else
    (void)threads;
    return 1;
#endif
}

// Captures the first exception thrown inside a parallel region so it can be
// rethrown once the region has closed.
class OmpError {
public:
    bool tripped() const { return failed_.load(std::memory_order_relaxed); }

    void capture() {
#pragma omp critical(hammock_omp_error)
        {
            if (!err_) err_ = std::current_exception();
        }
        failed_.store(true, std::memory_order_relaxed);
    }

    void rethrow_if_set() const {
        if (err_) std::rethrow_exception(err_);
    }

private:
    std::atomic<bool> failed_{false};
    std::exception_ptr err_;
};

// Compute pairwise Jaccard between two lists of HLL sketches into a (N, M) matrix.
py::array_t<double> pairwise_jaccard_hll(const std::vector<const HLLSketch*>& a,
                                         const std::vector<const HLLSketch*>& b,
                                         int threads) {
    require_uniform_precision(a, b);
    const int nt = pairwise_team_size(threads);
    const py::ssize_t n = static_cast<py::ssize_t>(a.size());
    const py::ssize_t m = static_cast<py::ssize_t>(b.size());
    py::array_t<double> out({n, m});
    auto buf = out.mutable_unchecked<2>();
    OmpError err;
    {
        py::gil_scoped_release release;
#pragma omp parallel for collapse(2) schedule(static) num_threads(nt)
        for (py::ssize_t i = 0; i < n; i++) {
            for (py::ssize_t j = 0; j < m; j++) {
                if (err.tripped()) continue;
                try {
                    buf(i, j) = a[i]->jaccard_similarity(*b[j]);
                } catch (...) {
                    err.capture();
                }
            }
        }
    }
    err.rethrow_if_set();
    return out;
}

// Compute pairwise Jaccard plus both directional containments in one pass.
// Returns (jaccard, containment_AB, containment_BA) where
//   containment_AB[i, j] = |a[i] ∩ b[j]| / |a[i]|
//   containment_BA[i, j] = |a[i] ∩ b[j]| / |b[j]|
// Co-sketch summaries (geom / arith / max) are derived Python-side from
// these two — keeps the binding narrow.
std::tuple<py::array_t<double>, py::array_t<double>, py::array_t<double>>
pairwise_metrics_hll(const std::vector<const HLLSketch*>& a,
                     const std::vector<const HLLSketch*>& b,
                     int threads) {
    require_uniform_precision(a, b);
    const int nt = pairwise_team_size(threads);
    const py::ssize_t n = static_cast<py::ssize_t>(a.size());
    const py::ssize_t m = static_cast<py::ssize_t>(b.size());
    py::array_t<double> jaccard({n, m});
    py::array_t<double> cont_ab({n, m});
    py::array_t<double> cont_ba({n, m});
    auto jbuf = jaccard.mutable_unchecked<2>();
    auto abbuf = cont_ab.mutable_unchecked<2>();
    auto babuf = cont_ba.mutable_unchecked<2>();
    OmpError err;
    {
        py::gil_scoped_release release;
        std::vector<double> a_card(n);
        std::vector<double> b_card(m);
#pragma omp parallel for schedule(static) num_threads(nt)
        for (py::ssize_t i = 0; i < n; i++) {
            a_card[i] = a[i]->cardinality();
        }
#pragma omp parallel for schedule(static) num_threads(nt)
        for (py::ssize_t j = 0; j < m; j++) {
            b_card[j] = b[j]->cardinality();
        }
#pragma omp parallel for collapse(2) schedule(static) num_threads(nt)
        for (py::ssize_t i = 0; i < n; i++) {
            for (py::ssize_t j = 0; j < m; j++) {
                if (err.tripped()) continue;
                try {
                    // One register pass yields both the Jaccard and the union
                    // cardinality. intersection_size() would re-estimate both
                    // operands' cardinalities per pair (they are already
                    // hoisted above) and materialize a union sketch — 16 MiB at
                    // p=24, allocated and freed N*M times. See
                    // HLLSketch::jaccard_and_union_cardinality for why this is
                    // bit-identical rather than merely close.
                    double u;
                    a[i]->jaccard_and_union_cardinality(*b[j], jbuf(i, j), u);
                    const double inter = std::max(0.0, a_card[i] + b_card[j] - u);
                    abbuf(i, j) = (a_card[i] > 0) ? (inter / a_card[i]) : 0.0;
                    babuf(i, j) = (b_card[j] > 0) ? (inter / b_card[j]) : 0.0;
                } catch (...) {
                    err.capture();
                }
            }
        }
    }
    err.rethrow_if_set();
    return {std::move(jaccard), std::move(cont_ab), std::move(cont_ba)};
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "hammock C++ core: HLL sketches and BED-mode processors.";

    py::class_<HLLSketch>(m, "HLLSketch")
        .def(py::init<size_t>(), py::arg("precision") = 18,
             "Construct an empty HyperLogLog sketch with the given precision (4..24).")
        .def("add_hash64",
             [](HLLSketch& self, uint64_t h) { self.add(h); },
             py::arg("hash_val"),
             "Add a precomputed 64-bit hash directly to the sketch.")
        .def("add_string",
             [](HLLSketch& self, const std::string& s) {
                 self.add(xxhash::hash64(s));
             },
             py::arg("s"),
             "Hash a string with XXH64(seed=0) and add to the sketch.")
        .def("estimate_cardinality", &HLLSketch::cardinality)
        .def("estimate_jaccard",
             [](const HLLSketch& self, const HLLSketch& other) {
                 return self.jaccard_similarity(other);
             },
             py::arg("other"))
        // Inclusion-exclusion |A| + |B| - |A u B|, clamped at 0 — the estimator
        // behind the containment/cosketch columns, and a different quantity
        // from estimate_jaccard's register-equality statistic. Exposed so the
        // two paths are separately testable; see
        // tests/test_containment_estimator.py.
        .def("estimate_intersection",
             [](const HLLSketch& self, const HLLSketch& other) {
                 return self.intersection_size(other);
             },
             py::arg("other"))
        // NOTE: a `merge_new` binding lived here until v0.6.0. Its only caller
        // was MinimizerSketch.merged(), which built the removed `_with_ends`
        // sketch (CLAUDE.md divergence #8). The underlying HLLSketch::union_with
        // is still used by pairwise_metrics_hll and hammock_cli.
        .def("clear", &HLLSketch::clear)
        .def("precision", &HLLSketch::precision)
        .def("__repr__", [](const HLLSketch& self) {
            return "<HLLSketch precision=" + std::to_string(self.precision()) + ">";
        });

    m.def("sketch_bed_file_hll", &sketch_bed_file_hll,
          py::arg("path"),
          py::arg("mode"),
          py::arg("precision") = 18,
          py::arg("separator") = "\t",
          py::arg("sub_a") = 1.0,
          py::arg("sub_b") = 1.0,
          py::arg("exp_a") = 0.0,
          py::arg("subB_method") = std::string("mixed-stride"),
          py::arg("seed") = 42,
          py::arg("gate_seed") = 31337,
          py::arg("verbose") = false,
          py::call_guard<py::gil_scoped_release>(),
          "Build an HLL sketch from one BED file in mode A, B, or C.");

    m.def("pairwise_jaccard_hll", &pairwise_jaccard_hll,
          py::arg("a").noconvert(), py::arg("b").noconvert(), py::arg("threads") = 0,
          "Compute the (len(a), len(b)) pairwise-Jaccard matrix for HLL sketches.\n"
          "threads=0 (default) leaves the OpenMP team size alone.");

    m.def("pairwise_metrics_hll", &pairwise_metrics_hll,
          py::arg("a").noconvert(), py::arg("b").noconvert(), py::arg("threads") = 0,
          "Compute (Jaccard, containment_AB, containment_BA) matrices in one pass.");

    m.def("read_filepath_list", &read_filepath_list, py::arg("list_file"));
}
