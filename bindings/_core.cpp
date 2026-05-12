#include "hammock/abstract_sketch.hpp"
#include "hammock/bagminhash_sketch.hpp"
#include "hammock/bed_parser.hpp"
#include "hammock/hll_sketch.hpp"
#include "hammock/processing_modes.hpp"
#include "hammock/xxhash.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <string>
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
                              int peak_height_column,
                              bool verbose) {
    const SubBMethod method = parse_subB_method(subB_method);
    HLLSketch sketch(precision);
    if (mode == "A") {
        process_bed_file_mode_a(path, sketch, seed, separator, peak_height_column, verbose);
    } else if (mode == "B") {
        process_bed_file_mode_b(path, sketch, seed, separator, sub_b, method,
                                gate_seed, peak_height_column, verbose);
    } else if (mode == "C") {
        process_bed_file_mode_c(path, sketch, seed, separator, sub_a, sub_b, exp_a,
                                method, gate_seed, peak_height_column, verbose);
    } else {
        throw std::invalid_argument("Mode must be one of A, B, C (got '" + mode + "')");
    }
    return sketch;
}

// Same, for a BagMinHash backing sketch.
BagMinHashSketch sketch_bed_file_bmh(const std::string& path,
                                     const std::string& mode,
                                     size_t num_hashes,
                                     uint64_t bmh_seed,
                                     const std::string& separator,
                                     double sub_a,
                                     double sub_b,
                                     double exp_a,
                                     const std::string& subB_method,
                                     uint64_t seed,
                                     uint32_t gate_seed,
                                     int peak_height_column,
                                     bool verbose) {
    const SubBMethod method = parse_subB_method(subB_method);
    BagMinHashSketch sketch(num_hashes, bmh_seed);
    if (mode == "A") {
        process_bed_file_mode_a(path, sketch, seed, separator, peak_height_column, verbose);
    } else if (mode == "B") {
        process_bed_file_mode_b(path, sketch, seed, separator, sub_b, method,
                                gate_seed, peak_height_column, verbose);
    } else if (mode == "C") {
        process_bed_file_mode_c(path, sketch, seed, separator, sub_a, sub_b, exp_a,
                                method, gate_seed, peak_height_column, verbose);
    } else {
        throw std::invalid_argument("Mode must be one of A, B, C (got '" + mode + "')");
    }
    return sketch;
}

// Compute pairwise Jaccard between two lists of HLL sketches into a (N, M) matrix.
py::array_t<double> pairwise_jaccard_hll(const std::vector<HLLSketch>& a,
                                         const std::vector<HLLSketch>& b) {
    const py::ssize_t n = static_cast<py::ssize_t>(a.size());
    const py::ssize_t m = static_cast<py::ssize_t>(b.size());
    py::array_t<double> out({n, m});
    auto buf = out.mutable_unchecked<2>();
    {
        py::gil_scoped_release release;
#pragma omp parallel for collapse(2) schedule(static)
        for (py::ssize_t i = 0; i < n; i++) {
            for (py::ssize_t j = 0; j < m; j++) {
                buf(i, j) = a[i].jaccard_similarity(b[j]);
            }
        }
    }
    return out;
}

// Compute pairwise Jaccard AND directional containment in one pass. Returns
// (jaccard, containment) where containment[i, j] = |a[i] ∩ b[j]| / |a[i]|.
// The expA exponent (mode C) is applied Python-side since it's CLI-metadata.
std::pair<py::array_t<double>, py::array_t<double>>
pairwise_jaccard_and_containment_hll(const std::vector<HLLSketch>& a,
                                     const std::vector<HLLSketch>& b) {
    const py::ssize_t n = static_cast<py::ssize_t>(a.size());
    const py::ssize_t m = static_cast<py::ssize_t>(b.size());
    py::array_t<double> jaccard({n, m});
    py::array_t<double> containment({n, m});
    auto jbuf = jaccard.mutable_unchecked<2>();
    auto cbuf = containment.mutable_unchecked<2>();
    {
        py::gil_scoped_release release;
        std::vector<double> self_card(n);
#pragma omp parallel for schedule(static)
        for (py::ssize_t i = 0; i < n; i++) {
            self_card[i] = a[i].cardinality();
        }
#pragma omp parallel for collapse(2) schedule(static)
        for (py::ssize_t i = 0; i < n; i++) {
            for (py::ssize_t j = 0; j < m; j++) {
                jbuf(i, j) = a[i].jaccard_similarity(b[j]);
                const double inter = a[i].intersection_size(b[j]);
                cbuf(i, j) = (self_card[i] > 0) ? (inter / self_card[i]) : 0.0;
            }
        }
    }
    return {std::move(jaccard), std::move(containment)};
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "hammock C++ core: HLL/BagMinHash sketches and BED-mode processors.";

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
        .def("merge_new",
             [](const HLLSketch& self, const HLLSketch& other) {
                 auto u = self.union_with(other);
                 // union_with returns unique_ptr<AbstractSketch>; we know the
                 // dynamic type is HLLSketch when both inputs are.
                 HLLSketch* hll_ptr = dynamic_cast<HLLSketch*>(u.get());
                 if (!hll_ptr) {
                     throw std::runtime_error("merge_new produced a non-HLL sketch");
                 }
                 return HLLSketch(*hll_ptr);
             },
             py::arg("other"),
             "Return a new HLLSketch that is the union of self and other.")
        .def("clear", &HLLSketch::clear)
        .def("precision", &HLLSketch::precision)
        .def("__repr__", [](const HLLSketch& self) {
            return "<HLLSketch precision=" + std::to_string(self.precision()) + ">";
        });

    py::class_<BagMinHashSketch>(m, "BagMinHashSketch")
        .def(py::init<size_t, uint64_t>(),
             py::arg("num_hashes"), py::arg("seed") = 0)
        .def("add_hash64",
             [](BagMinHashSketch& self, uint64_t h) { self.add(h); },
             py::arg("hash_val"))
        .def("add_with_normalized_count",
             &BagMinHashSketch::add_with_normalized_count,
             py::arg("hash_val"), py::arg("raw_count"), py::arg("scale_factor") = 100.0)
        .def("estimate_cardinality", &BagMinHashSketch::cardinality)
        .def("estimate_jaccard",
             [](const BagMinHashSketch& self, const BagMinHashSketch& other) {
                 return self.jaccard_similarity(other);
             },
             py::arg("other"))
        .def("clear", &BagMinHashSketch::clear);

    m.def("sketch_bed_file_hll", &sketch_bed_file_hll,
          py::arg("path"),
          py::arg("mode"),
          py::arg("precision") = 18,
          py::arg("separator") = "\t",
          py::arg("sub_a") = 1.0,
          py::arg("sub_b") = 1.0,
          py::arg("exp_a") = 0.0,
          py::arg("subB_method") = std::string("hash-threshold"),
          py::arg("seed") = 42,
          py::arg("gate_seed") = 31337,
          py::arg("peak_height_column") = -1,
          py::arg("verbose") = false,
          py::call_guard<py::gil_scoped_release>(),
          "Build an HLL sketch from one BED file in mode A, B, or C.");

    m.def("sketch_bed_file_bmh", &sketch_bed_file_bmh,
          py::arg("path"),
          py::arg("mode"),
          py::arg("num_hashes") = 128,
          py::arg("bmh_seed") = 0,
          py::arg("separator") = "\t",
          py::arg("sub_a") = 1.0,
          py::arg("sub_b") = 1.0,
          py::arg("exp_a") = 0.0,
          py::arg("subB_method") = std::string("hash-threshold"),
          py::arg("seed") = 42,
          py::arg("gate_seed") = 31337,
          py::arg("peak_height_column") = -1,
          py::arg("verbose") = false,
          py::call_guard<py::gil_scoped_release>(),
          "Build a BagMinHash sketch from one BED file (use this when --count <col> is set).");

    m.def("pairwise_jaccard_hll", &pairwise_jaccard_hll,
          py::arg("a"), py::arg("b"),
          "Compute the (len(a), len(b)) pairwise-Jaccard matrix for HLL sketches.");

    m.def("pairwise_jaccard_and_containment_hll", &pairwise_jaccard_and_containment_hll,
          py::arg("a"), py::arg("b"),
          "Compute (Jaccard, directional containment) matrices in one pass.");

    m.def("read_filepath_list", &read_filepath_list, py::arg("list_file"));
}
