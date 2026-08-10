#!/bin/bash
#SBATCH --job-name=fig3_lowN
#SBATCH --partition=parallel
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH --output=experiments/bedtools_benchmark/logs/fig3_lowN_%j.out
#SBATCH --error=experiments/bedtools_benchmark/logs/fig3_lowN_%j.err

# Panel A's low-N cells (N=2..32) are unreliable at the main run's n=3: a local,
# unallocated, paired 20/15-rep check (same code path, same session, not this
# job) found bedtools wins 20/20 at N=2 (median 0.167s vs hammock 0.256s) and
# N=2->N=4 is monotonic (0.167 -> 0.486s) -- both contradicting the published
# n=3 CSV, which read N=2 hammock-faster (2.02x) and a physically-impossible
# N=2->N=4 cost DROP (0.50 -> 0.18s). That check was good enough to prove the
# published low-N cells are wrong, but was run on a shared/interactive node,
# not this repo's usual exclusive allocation -- not the number to commit.
#
# This job re-measures N in {2,4,8,16,32} at n=20 on its own exclusive
# allocation, so the replacement rows are measured the same way the rest of
# Panel A was. N>=64 rows are untouched -- their CV was already low (see
# docs/data/cpp_vs_bedtools_t16_p18.csv) and 20 reps there would cost far more
# for no benefit (N=512 alone would balloon to ~2.4 hours at n=20).

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude
ml bedtools/2.30.0
ml parallel
mkdir -p experiments/bedtools_benchmark/logs

PYTHON=/home/jbonnie1/.conda/envs/hammock/bin/python3
export HAMMOCK_CPP_BIN=/home/jbonnie1/.conda/envs/claude-ref-comparison/lib/python3.10/site-packages/bin/hammock-cpp

echo "node: $(hostname)  cores: $(nproc)"
"$HAMMOCK_CPP_BIN" --version

"$PYTHON" experiments/bedtools_benchmark/benchmark_cpp_vs_bedtools.py \
    --threads 16 \
    --precision 18 \
    --num-files 2,4,8,16,32 \
    --runs 20 \
    --corpus-seed 20260810 \
    "$@"
