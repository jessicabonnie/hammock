#ifndef HAMMOCK_OMP_UTIL_HPP
#define HAMMOCK_OMP_UTIL_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

// Team size for an OpenMP region driven by the CLI's `--threads`/`threads=`
// value. `threads <= 0` means "whatever OpenMP would have picked", so passing
// omp_get_max_threads() back in is a no-op -- this is what `hammock_cli.cpp`
// relies on: it sets the ambient team size once, globally, via
// omp_set_num_threads() before its per-file sketching loop and its pairwise
// loop, and passes threads=0 into everything below it so those regions just
// inherit that ambient value unchanged.
//
// Originally existed only for the pairwise loops (bindings/_core.cpp): without
// it the pairwise phase ignored --threads entirely, since omp_set_num_threads
// appears only in hammock_cli.cpp, so `hammock --threads 4` inside a 4-CPU
// cgroup still ran that loop on every core on the node. Extended to the
// sketching-phase regions in processing_modes.cpp to close the matching gap
// there -- see docs/seed-hammock-cpp-file-dispatch.md Part 1. Previously each
// of the two call sites carried its own hand-copied version of this function
// (`pairwise_team_size` / `sketch_team_size`); hoisted here so the clamp
// formula can't silently drift between the two front ends' translation units.
//
// Applied as a num_threads() clause rather than omp_set_num_threads() in both
// cases so it cannot leak into unrelated regions.
inline int omp_team_size(int threads) {
#ifdef _OPENMP
    // Clamp rather than trust: --threads is an unvalidated `type=int` on the
    // CLI, and since it began driving num_threads() a fat-fingered --threads
    // 100000 would ask libgomp for 100000 threads instead of merely oversizing
    // a Python pool. 4x the machine's max is well past any real use.
    const int cap = 4 * omp_get_num_procs();
    if (threads > cap) return cap;
    return (threads > 0) ? threads : omp_get_max_threads();
#else
    (void)threads;
    return 1;
#endif
}

#endif  // HAMMOCK_OMP_UTIL_HPP
