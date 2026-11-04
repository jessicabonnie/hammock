#!/usr/bin/env bash
#SBATCH --job-name=seq-repro-20261028
#SBATCH --partition=parallel
#SBATCH --account=blangme2_scal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --exclude=bigmem23,c189,c712,c717
#SBATCH --output=/vast/blangme2/jbonnie/hammock-seq-repro-20261028/results/logs/slurm/%A_%a.out
#SBATCH --error=/vast/blangme2/jbonnie/hammock-seq-repro-20261028/results/logs/slurm/%A_%a.err
# Reproducibility check: rerun archived Mode D (rehash-selector64) jobs with an
# isolated build of current main and write ONLY under the new result root.
set -euo pipefail

CHECK_DIR=/home/jbonnie1/interval_sketch/hammock_claude/experiments/seq_repro_check_20261028
BUILD=/vast/blangme2/jbonnie/hammock-seq-repro-20261028
RESULT_ROOT="$BUILD/results"
ENV=/home/jbonnie1/.conda/envs/claude-ref-comparison
ARCHIVED_RUNTIME=/vast/blangme2/jbonnie/hammock-mode-d-rehash-sweep.sCydsC/runtime/site-packages

task=${SLURM_ARRAY_TASK_ID:?array task required}
line=$(awk -F'\t' -v t="$task" 'NR>1 && $1==t' "$CHECK_DIR/jobs.tsv")
[[ -n $line ]] || { echo "no job row for task $task" >&2; exit 2; }
IFS=$'\t' read -r _ label runtime precision seed k w fasta_list archived_csv <<< "$line"

case $runtime in
  new) export PYTHONPATH="$BUILD/runtime/site-packages" ;;
  archived) export PYTHONPATH="$ARCHIVED_RUNTIME" ;;
  *) echo "bad runtime $runtime" >&2; exit 2 ;;
esac
export LD_PRELOAD="$ENV/lib/libstdc++.so.6"

extra=()
[[ $label == fig5_* ]] && extra+=(--full-paths)

out_dir="$RESULT_ROOT/task$(printf %02d "$task")_${label}_${runtime}_p${precision}_seed${seed}_k${k}_w${w}"
if [[ -e $out_dir ]]; then
  echo "refusing to reuse existing output dir: $out_dir" >&2; exit 3
fi
mkdir -p "$out_dir"
cd "$out_dir"

"$ENV/bin/python" -c 'import hammock, hammock._core, digest, xxhash, sys; print("hammock", hammock.__file__); print("_core", hammock._core.__file__); print("digest", digest.__file__); print("xxhash", xxhash.VERSION); print("python", sys.version.split()[0])' > provenance.txt

/usr/bin/time -v -o time.txt "$ENV/bin/python" -m hammock "$fasta_list" "$fasta_list" \
  --mode D --sequence-hll-hash rehash-selector64 \
  --seed "$seed" --precision "$precision" -k "$k" -w "$w" \
  --threads 1 --metrics "${extra[@]}" \
  --outprefix "$out_dir/hammock" > stdout.txt 2> stderr.txt

ls "$out_dir"/*.csv
echo "done task=$task label=$label runtime=$runtime"
