#!/bin/bash
#SBATCH --job-name=maurano_dhs_dag
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=experiments/maurano_dhs_validation/logs/slurm/dag_%j.out
#SBATCH --error=experiments/maurano_dhs_validation/logs/slurm/dag_%j.err

# Drives the Snakemake DAG for maurano_dhs_validation. This job stays
# resident as the controller; it submits per-rule sbatch jobs via the
# workflow/slurm_profile/. Heavy compute happens in child jobs, so this
# parent only needs a single cpu and a long time budget.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude/experiments/maurano_dhs_validation

CONDA_ENV=/home/jbonnie1/.conda/envs/claude-ref-comparison
SNAKEMAKE="$CONDA_ENV/bin/snakemake"

mkdir -p logs/slurm

# Sanity: print the planned DAG into the log before launching.
"$SNAKEMAKE" -s workflow/Snakefile -n --quiet | tail -n 20

# Launch the workflow. --profile picks up workflow/slurm_profile/config.yaml.
"$SNAKEMAKE" -s workflow/Snakefile --profile workflow/slurm_profile/
