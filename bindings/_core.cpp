#include "hammock/abstract_sketch.hpp"
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/omp_util.hpp"
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
//
// `threads` bounds the OMP team inside process_bed_file_mode_b/c's sketching
// region (mode A has no parallel region) -- passed through unresolved, same
// convention hammock_cli.cpp uses for its own calls: process_bed_file_mode_b/c
// resolve it via omp_team_size() themselves exactly once, so this function
// deliberately does *not* pre-resolve it here too (an earlier version did,
// which meant a raw request got clamped twice through two copies of the same
// formula -- harmless while they agreed, but a latent drift hazard). Callers
// dispatching multiple files concurrently via their own thread pool (runner.py's
// `_sketch_many`) are expected to divide their overall thread budget across
// that pool's actual concurrency before calling in here, so the *total*
// concurrent OS threads stays bounded by the user's `--threads` -- see
// docs/seed-hammock-cpp-file-dispatch.md Part 1. threads<=0 leaves the team
// size at OpenMP's ambient default, same convention as the pairwise loops.
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
                              bool verbose,
                              int threads) {
    const SubBMethod method = parse_subB_method(subB_method);
    HLLSketch sketch(precision);
    if (mode == "A") {
        process_bed_file_mode_a(path, sketch, seed, separator, verbose);
    } else if (mode == "B") {
        process_bed_file_mode_b(path, sketch, seed, separator, sub_b, method,
                                gate_seed, verbose, threads);
    } else if (mode == "C") {
        process_bed_file_mode_c(path, sketch, seed, separator, sub_a, sub_b, exp_a,
                                method, gate_seed, verbose, threads);
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

// Borrow the sketches for the call, but hold a strong reference to each element
// while doing it.
//
// Taking `std::vector<const HLLSketch*>` via pybind's list_caster avoids the
// deep copy, but it also drops every guarantee that kept the by-value form
// safe. The args tuple keeps the *list object* alive, not its elements, and
// list_caster's element accessor releases each item as it iterates. Two
// consequences, both reproduced as hard crashes before this helper existed:
//
//   * another Python thread calling lst.clear() (or del lst[:], lst[i] = x)
//     during the GIL-released region drops the last reference to every sketch
//     and frees it mid-loop -- SIGSEGV, no traceback;
//   * any non-list sequence whose __getitem__ synthesizes objects (the "slow
//     sequence" path list_caster accepts) has no other owner at all, so every
//     pointer dangles immediately -- SIGABRT.
//
// One incref per element restores the old safety at O(N) refcount operations
// instead of O(N * 2^p) bytes copied, which is the whole point of borrowing.
static void borrow_sketches(const py::sequence& seq,
                            const char* argname,
                            std::vector<py::object>& keepalive,
                            std::vector<const HLLSketch*>& out) {
    const py::ssize_t n = static_cast<py::ssize_t>(py::len(seq));
    keepalive.reserve(static_cast<size_t>(n));
    out.reserve(static_cast<size_t>(n));
    for (py::ssize_t i = 0; i < n; ++i) {
        py::object item = seq[static_cast<size_t>(i)];
        if (item.is_none()) {
            throw py::type_error(std::string(argname) + "[" + std::to_string(i) +
                                 "] is None, expected an HLLSketch");
        }
        const HLLSketch* ptr = nullptr;
        try {
            ptr = item.cast<const HLLSketch*>();
        } catch (const py::cast_error&) {
            ptr = nullptr;
        }
        if (!ptr) {
            throw py::type_error(std::string(argname) + "[" + std::to_string(i) +
                                 "] is not an HLLSketch");
        }
        keepalive.push_back(std::move(item));
        out.push_back(ptr);
    }
}

// Compute pairwise Jaccard between two lists of HLL sketches into a (N, M) matrix.
py::array_t<double> pairwise_jaccard_hll(const py::sequence& a_seq,
                                         const py::sequence& b_seq,
                                         int threads) {
    std::vector<py::object> a_keep, b_keep;   // outlive the GIL-released region
    std::vector<const HLLSketch*> a, b;
    borrow_sketches(a_seq, "a", a_keep, a);
    borrow_sketches(b_seq, "b", b_keep, b);
    require_uniform_precision(a, b);
    const int nt = omp_team_size(threads);
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
                    buf(i, j) = a[i]->reg_eq_similarity(*b[j]);
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
pairwise_metrics_hll(const py::sequence& a_seq,
                     const py::sequence& b_seq,
                     int threads) {
    std::vector<py::object> a_keep, b_keep;   // outlive the GIL-released region
    std::vector<const HLLSketch*> a, b;
    borrow_sketches(a_seq, "a", a_keep, a);
    borrow_sketches(b_seq, "b", b_keep, b);
    require_uniform_precision(a, b);
    const int nt = omp_team_size(threads);
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
        // Latched like the pair loop below, not bare. cardinality() allocates
        // its histogram, so std::bad_alloc is reachable here, and an exception
        // escaping an OpenMP structured block calls std::terminate -- the
        // interpreter dies with no traceback. The comment on
        // require_uniform_precision claims OmpError catches "whatever else might
        // still throw"; before this it did not cover these two loops.
#pragma omp parallel for schedule(static) num_threads(nt)
        for (py::ssize_t i = 0; i < n; i++) {
            if (err.tripped()) continue;
            try {
                a_card[i] = a[i]->cardinality();
            } catch (...) {
                err.capture();
            }
        }
#pragma omp parallel for schedule(static) num_threads(nt)
        for (py::ssize_t j = 0; j < m; j++) {
            if (err.tripped()) continue;
            try {
                b_card[j] = b[j]->cardinality();
            } catch (...) {
                err.capture();
            }
        }
#pragma omp parallel for collapse(2) schedule(static) num_threads(nt)
        for (py::ssize_t i = 0; i < n; i++) {
            for (py::ssize_t j = 0; j < m; j++) {
                if (err.tripped()) continue;
                try {
                    // One register pass yields both the register-equality
                    // similarity and the union cardinality. intersection_size()
                    // would re-estimate both operands' cardinalities per pair
                    // (they are already hoisted above) and materialize a union
                    // sketch — 16 MiB at p=24, allocated and freed N*M times. See
                    // HLLSketch::reg_eq_and_union_cardinality for why this is
                    // bit-identical rather than merely close.
                    double u;
                    a[i]->reg_eq_and_union_cardinality(*b[j], jbuf(i, j), u);
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
        .def("estimate_reg_eq_similarity",
             [](const HLLSketch& self, const HLLSketch& other) {
                 return self.reg_eq_similarity(other);
             },
             py::arg("other"))
        // Inclusion-exclusion |A| + |B| - |A u B|, clamped at 0 — the estimator
        // behind the containment/cosketch columns, and a different quantity
        // from estimate_reg_eq_similarity's register-equality statistic. Exposed so the
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
          py::arg("threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Build an HLL sketch from one BED file in mode A, B, or C.\n"
          "threads bounds the sketching-phase OMP team (modes B/C only); "
          "threads<=0 (default) leaves it at OpenMP's ambient default. "
          "Callers dispatching many files concurrently via their own thread "
          "pool should divide their overall thread budget by that pool's "
          "actual concurrency before passing it in here.");

    m.def("pairwise_jaccard_hll", &pairwise_jaccard_hll,
          py::arg("a"), py::arg("b"), py::arg("threads") = 0,
          "Compute the (len(a), len(b)) pairwise-Jaccard matrix for HLL sketches.\n"
          "threads=0 (default) leaves the OpenMP team size alone.");

    m.def("pairwise_metrics_hll", &pairwise_metrics_hll,
          py::arg("a"), py::arg("b"), py::arg("threads") = 0,
          "Compute (Jaccard, containment_AB, containment_BA) matrices in one pass.");

    m.def("read_filepath_list", &read_filepath_list, py::arg("list_file"));
}
