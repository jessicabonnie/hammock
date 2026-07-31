#!/usr/bin/env python3
"""Driver: run the R column-comparison script via subprocess.

R's hclust complete-linkage and scipy's break ties differently on these
tightly-bunched similarities; the existing pipeline (cluster_plot.R) and
its headline ARIs use R's hclust, so we use R to reproduce them.
"""
import subprocess
import sys

RSCRIPT = "/data/apps/extern/spack_on/gcc/9.3.0/r/4.3.0-5w26bclf7ogtbbwwfx337icyt6p53t3m/bin/Rscript"
R_SCRIPT = "/home/jbonnie1/interval_sketch/hammock_claude/experiments/mus-homo/scripts/compute_column_comparison.R"

p = subprocess.run([RSCRIPT, R_SCRIPT], capture_output=True, text=True)
sys.stdout.write(p.stdout)
sys.stderr.write(p.stderr)
sys.exit(p.returncode)
