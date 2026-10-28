#!/bin/bash
#SBATCH --job-name=hm_build_6248r
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=experiments/vectorized_rebuild_20261027/logs/build_%j.out
#SBATCH --error=experiments/vectorized_rebuild_20261027/logs/build_%j.err

# Build the db1411a hammock-cpp CLI on a Xeon Gold 6248R benchmark node with
# the conda-forge GCC 14.3 toolchain. The Spack GCC 9.3 module (on PATH by
# default) produced a largely unvectorized binary that sketched 27--44% slower.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude

EXP=experiments/vectorized_rebuild_20261027
WT=/vast/blangme2/jbonnie/hammock-db1411a-c6248r
BUILD=/vast/blangme2/jbonnie/hammock-db1411a-builds/c6248r-gcc14
ENVBIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin
CMAKE=$ENVBIN/cmake
CXX=$ENVBIN/x86_64-conda-linux-gnu-g++
CC=$ENVBIN/x86_64-conda-linux-gnu-gcc

cpu=$(lscpu | sed -n 's/^Model name:[[:space:]]*//p')
[[ "$cpu" == *"6248R"* ]] || { echo "unexpected CPU: $cpu" >&2; exit 1; }
[[ "$(git -C "$WT" rev-parse HEAD)" == db1411a* ]] || { echo "worktree not at db1411a" >&2; exit 1; }
[[ -z "$(git -C "$WT" status --porcelain)" ]] || { echo "worktree dirty" >&2; exit 1; }

"$CMAKE" -S "$WT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" \
  -DSKBUILD_PROJECT_VERSION=0.12.0
"$CMAKE" --build "$BUILD" --target hammock-cpp -j 16

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
} | tee "$EXP/build_metadata.txt"

(( ymm > 100 )) || { echo "binary is not vectorized (ymm=$ymm)" >&2; exit 1; }
