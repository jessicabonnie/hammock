#!/usr/bin/env bash
#SBATCH --job-name=xref-rehash-p24
#SBATCH --partition=parallel
#SBATCH --account=blangme2_scal
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=bigmem23,c189,c712,c717
#SBATCH --output=/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/logs/slurm/%A_%a.out
#SBATCH --error=/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64/logs/slurm/%A_%a.err

set -euo pipefail

REPO=/home/jbonnie1/interval_sketch/hammock_claude
EXPERIMENT="$REPO/experiments/ref-comparison/rehash_rerun"
RESULT_ROOT=/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64
FASTA_ROOT=/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/fastas
HAMMOCK=/home/jbonnie1/.conda/envs/claude-ref-comparison/bin/hammock
SEEDS=(1 42 99)
PEAK_TYPES=(broad narrow)
CELL_COUNT=63
TASKS_PER_SEED=$((CELL_COUNT * ${#PEAK_TYPES[@]}))

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
if (( task_id < 0 || task_id >= ${#SEEDS[@]} * TASKS_PER_SEED )); then
    echo "array task $task_id is outside the expected range 0-377" >&2
    exit 2
fi

seed_index=$((task_id / TASKS_PER_SEED))
within_seed=$((task_id % TASKS_PER_SEED))
peak_index=$((within_seed / CELL_COUNT))
cell_index=$((within_seed % CELL_COUNT + 2))
seed=${SEEDS[$seed_index]}
peak_type=${PEAK_TYPES[$peak_index]}
read -r k w < <(sed -n "${cell_index}p" "$EXPERIMENT/cells.tsv")

if [[ -z ${k:-} || -z ${w:-} ]]; then
    echo "failed to resolve cell index $cell_index" >&2
    exit 2
fi

cell_dir="$RESULT_ROOT/seed${seed}/${peak_type}/k${k}_w${w}"
log_dir="$RESULT_ROOT/logs/seed${seed}/${peak_type}"
mkdir -p "$cell_dir" "$log_dir"
fasta_list="$cell_dir/exp_a_all_fastas.txt"
: > "$fasta_list"
for ref in CHM13 GRCh37 GRCh38; do
    for sample in ENCSR175ABH_H3K27ac ENCSR864OOO_H3K27ac ENCSR954JMZ_H3K27ac; do
        fasta="$FASTA_ROOT/$peak_type/$ref/$sample.fa"
        [[ -s $fasta ]] || { echo "missing FASTA: $fasta" >&2; exit 3; }
        printf '%s\n' "$fasta" >> "$fasta_list"
    done
done

output="$cell_dir/exp_a_mnmzr_p24_jaccD_k${k}_w${w}_rehash-selector64_full.csv"
if [[ -s $output ]]; then
    rows=$(wc -l < "$output")
    if [[ $rows -eq 82 ]]; then
        echo "validated existing output: $output"
        exit 0
    fi
    echo "refusing incomplete existing output ($rows lines): $output" >&2
    exit 4
fi

export PYTHONPATH="$REPO/python${PYTHONPATH:+:$PYTHONPATH}"
"$HAMMOCK" "$fasta_list" "$fasta_list" \
    --mode D \
    --outprefix "$cell_dir/exp_a" \
    --sequence-hll-hash rehash-selector64 \
    --seed "$seed" \
    --precision 24 \
    -k "$k" \
    -w "$w" \
    --threads 1 \
    --full-paths \
    --metrics \
    2> "$log_dir/k${k}_w${w}.log"

rows=$(wc -l < "$output")
[[ $rows -eq 82 ]] || { echo "expected 82 CSV lines, found $rows: $output" >&2; exit 5; }
printf 'completed seed=%s peak_type=%s k=%s w=%s output=%s\n' \
    "$seed" "$peak_type" "$k" "$w" "$output"
