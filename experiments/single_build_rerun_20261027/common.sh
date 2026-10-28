# Sourced by the sbatch_*.sh wrappers in this directory. Not run directly.
#
# Single-build rerun (2026-10-27): every hammock timing in this directory is
# taken with ONE hammock-cpp 0.12.0 binary (commit db1411a, conda-forge GCC 14.3,
# -O3 -march=native, built on a Xeon Gold 6248R node by job 32614423), except
# the per-pair cost rerun on the 6448Y node type, which needs a native rebuild
# of the same commit with the same compiler (sbatch_build_6448y.sh).
#
# Nothing here writes outside $EXP.

EXP=/home/jbonnie1/interval_sketch/hammock_claude/experiments/single_build_rerun_20261027
REPO=/home/jbonnie1/interval_sketch/hammock_claude
SINGLE_BINARY=/vast/blangme2/jbonnie/hammock-db1411a-builds/c6248r-gcc14/hammock-cpp
SINGLE_SHA=c6c13e8e469404cac23b538a84f71b4fad04ccaf628dcd76ec1c0f1a836d3953
PY_HAMMOCK=/home/jbonnie1/.conda/envs/hammock/bin/python3
PY_REF=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/python3

require_cpu() {
  local want=$1 cpu
  cpu=$(lscpu | sed -n 's/^Model name:[[:space:]]*//p')
  [[ "$cpu" == *"$want"* ]] || { echo "unexpected CPU: $cpu (want $want)" >&2; exit 1; }
}

require_sha() {
  local bin=$1 want=$2 got
  got=$(sha256sum "$bin" | cut -d' ' -f1)
  [[ "$got" == "$want" ]] || { echo "binary sha256 mismatch for $bin: $got" >&2; exit 1; }
}

# write_metadata <file> <binary> ; appends nothing if the file exists (refuse).
write_metadata() {
  local out=$1 bin=$2
  [[ ! -e "$out" ]] || { echo "refusing to overwrite $out" >&2; exit 1; }
  {
    echo "job_id: ${SLURM_JOB_ID:-none}"
    echo "job_name: ${SLURM_JOB_NAME:-none}"
    echo "partition: ${SLURM_JOB_PARTITION:-none}"
    echo "node: $(hostname)"
    echo "cpu: $(lscpu | sed -n 's/^Model name:[[:space:]]*//p')"
    echo "cores_visible: $(nproc)"
    echo "slurm_cpus_on_node: ${SLURM_CPUS_ON_NODE:-unset}"
    echo "slurm_cpus_per_task: ${SLURM_CPUS_PER_TASK:-unset}"
    echo "tmpdir: ${TMPDIR:-/tmp} ($(df -hT "${TMPDIR:-/tmp}" | tail -1 | awk '{print $2, $3, $5}'))"
    echo "harness_repo_commit: $(git -C "$REPO" rev-parse HEAD)$( [[ -n "$(git -C "$REPO" status --porcelain)" ]] && echo ' (repo working tree dirty; harness copies in harness/)')"
    echo "harness_copy_sha256:"
    (cd "$EXP/harness" && sha256sum benchmark_cpp_vs_bedtools.py sweep.py pairwise_cost_by_precision.py bedtools.sh | sed 's/^/  /')
    echo "binary: $bin"
    echo "binary_sha256: $(sha256sum "$bin" | cut -d' ' -f1)"
    echo "binary_version: $("$bin" --version)"
    echo "binary_libs:"
    ldd "$bin" | sed 's/^/  /'
    echo "bedtools_version: $(command -v bedtools >/dev/null && bedtools --version || echo 'not loaded')"
    echo "bedtools_path: $(command -v bedtools || echo none)"
    echo "parallel_version: $(command -v parallel >/dev/null && parallel --version | head -1 || echo 'not loaded')"
    echo "python: $("$PY_HAMMOCK" --version 2>&1) ($PY_HAMMOCK); $("$PY_REF" --version 2>&1) ($PY_REF)"
    echo "start: $(date -Is)"
  } | tee "$out"
}

finish_metadata() {
  echo "end: $(date -Is)" | tee -a "$1"
}
