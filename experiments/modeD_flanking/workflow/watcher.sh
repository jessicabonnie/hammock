#!/usr/bin/env bash
#SBATCH --job-name=mdf_watcher
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking/logs/watcher.%j.log

set -u
EXPT=/home/jbonnie1/interval_sketch/hammock_claude/experiments/modeD_flanking
cd "$EXPT"
echo "[watcher] started $(date)"

# Wait for both Maurano DAG (snakejob.maurano + snakejob.run_mode_d) and
# part2_sweep to drain.
while :; do
  pending=$(squeue -h -u "$USER" -o '%j' 2>/dev/null \
            | grep -cE '^(snakejob\.run_mode_d|snakejob\.maurano|maurano_dhs_dag|part2_sweep|snakejob\.analyze)' || true)
  if [ "${pending:-0}" = "0" ]; then
    break
  fi
  echo "[watcher] $(date) — $pending relevant jobs still queued/running; sleep 120"
  sleep 120
done

echo "[watcher] $(date) — both reruns complete; running Part 1 + Part 2 R analyses"
. /etc/profile 2>/dev/null || true
ml gcc/9.3.0 r/4.3.0 libjpeg/9c

echo "=== Part 1 ==="
Rscript analyze_part1_maurano.R 2>&1 | tail -8

echo "=== Part 2 ==="
Rscript analyze_part2_synthetic.R 2>&1 | tail -8

echo "[watcher] $(date) — done"
