#!/bin/bash
#SBATCH --job-name=sbr_build_6448y
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/build_6448y_%j.out
#SBATCH --error=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027/logs/build_6448y_%j.err

# Native db1411a build for the Xeon Gold 6448Y node type (Rockfish `shared`
# partition), needed only for the per-pair estimator-cost rerun, which the
# archived 0.7.0 measurement (job 29620391) took on a 6448Y node (sr09) with 16
# reserved cores. Copy of experiments/vectorized_rebuild_20261027/sbatch_build.sh
# with: CPU check 6448Y, shared partition with 16 reserved cores, and a NEW
# build directory inside this experiment directory. Same worktree (db1411a,
# clean), same conda-forge GCC 14.3 toolchain, same CMake options.

set -euo pipefail
EXP=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027
WT=/vast/blangme2/jbonnie/hammock-db1411a-c6248r
BUILD=$EXP/build_sr6448y_gcc14
META=$EXP/build_metadata_sr6448y.txt
ENVBIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin
CMAKE=$ENVBIN/cmake
CXX=$ENVBIN/x86_64-conda-linux-gnu-g++
CC=$ENVBIN/x86_64-conda-linux-gnu-gcc

[[ ! -e "$BUILD" ]] || { echo "build dir $BUILD already exists; refusing to reuse" >&2; exit 1; }
[[ ! -e "$META" ]] || { echo "$META already exists; refusing to overwrite" >&2; exit 1; }
cpu=$(lscpu | sed -n 's/^Model name:[[:space:]]*//p')
[[ "$cpu" == *"6448Y"* ]] || { echo "unexpected CPU: $cpu" >&2; exit 1; }
[[ "$(git -C "$WT" rev-parse HEAD)" == db1411a* ]] || { echo "worktree not at db1411a" >&2; exit 1; }
[[ -z "$(git -C "$WT" status --porcelain)" ]] || { echo "worktree dirty" >&2; exit 1; }

"$CMAKE" -S "$WT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" \
  -DSKBUILD_PROJECT_VERSION=0.12.0
"$CMAKE" --build "$BUILD" --target hammock-cpp -j 16

[[ -z "$(git -C "$WT" status --porcelain)" ]] || { echo "worktree became dirty during build" >&2; exit 1; }

BIN=$BUILD/hammock-cpp
zmm=$(objdump -d "$BIN" | grep -c '%zmm' || true)
ymm=$(objdump -d "$BIN" | grep -c '%ymm' || true)
{
  echo "job_id: ${SLURM_JOB_ID}"
  echo "build_node: $(hostname)"
  echo "build_cpu: $cpu"
  echo "commit: $(git -C "$WT" rev-parse HEAD)"
  echo "compiler: $("$CXX" --version | head -1)"
  echo "binary_comment: $(readelf -p .comment "$BIN" | grep -o 'GCC: (conda-forge[^)]*) [0-9.]*' | head -1)"
  echo "cxx_flags: $(grep -h CXX_FLAGS "$BUILD"/CMakeFiles/hammock-cpp.dir/flags.make)"
  echo "binary: $BIN"
  echo "binary_sha256: $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "binary_version: $("$BIN" --version)"
  echo "zmm_instructions: $zmm"
  echo "ymm_instructions: $ymm"
} | tee "$META"

(( ymm > 100 )) || { echo "binary is not vectorized (ymm=$ymm)" >&2; exit 1; }
