#!/bin/bash
#SBATCH --job-name=mus_homo_dag
#SBATCH --partition=shared
#SBATCH --account=blangme2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=experiments/mus-homo/logs/slurm/dag_%j.out
#SBATCH --error=experiments/mus-homo/logs/slurm/dag_%j.err

# Drives the Snakemake DAG for mus-homo. Modeled on
# experiments/maurano_dhs_validation/sbatch_workflow.sh.

set -euo pipefail
cd /home/jbonnie1/interval_sketch/hammock_claude/experiments/mus-homo

CONDA_ENV=/home/jbonnie1/.conda/envs/claude-ref-comparison
SNAKEMAKE="$CONDA_ENV/bin/snakemake"

mkdir -p logs/slurm

"$SNAKEMAKE" -s workflow/Snakefile -n --quiet | tail -n 20
"$SNAKEMAKE" -s workflow/Snakefile --profile workflow/slurm_profile/
