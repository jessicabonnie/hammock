# Results

The corrected cross-reference rerun completed on 2026-08-25. All 378 matrices
(63 cells x 2 peak types x 3 seeds) passed the analyzer's row-count, unique-pair,
symmetry, unit-diagonal, mode, precision, parameter, finiteness, and bounds
checks.

## Outcome

`k=20,w=20` remained the unique inclusion--exclusion separation leader for
both broad and narrow peaks at every prespecified seed.

| Peak type | Seed-42 delta | Three-seed range | Best rank | Worst rank |
|---|---:|---:|---:|---:|
| broad | 0.559013 | 0.558953--0.559156 | 1 | 1 |
| narrow | 0.595431 | 0.595214--0.595431 | 1 | 1 |

The largest three-seed delta range for any tested cell was 0.000510
(`broad`, `k=35,w=35`). At seed 42, the next broad cells were `k20_w25`
(0.557930) and `k20_w30` (0.557104); the lead is small but reproduced at all
three seeds. The nine seed-42 broad peak sets retain three monophyletic tissue
clades, each containing GRCh37, GRCh38, and CHM13.

The legacy-to-corrected seed-42 change at `k20_w20` was 0.563200 to 0.559013
for broad peaks and 0.598568 to 0.595431 for narrow peaks. The scientific
conclusion and selected cell survive; the paper-facing numerical value changes.

## Provenance

- Source implementation at launch: `18002b1c4bc2ecdfdb080f78cee2b6680228bdb8`.
- A concurrent commit changed paper figure sources during execution but did
  not change `python/hammock/modes/sequence.py`, `runner.py`, or `cli.py`.
- Hammock: 0.11.0; Python: the `claude-ref-comparison` environment.
- Parameters: Mode D, p=24, `rehash-selector64`, seeds 1/42/99, one thread.
- Slurm gate job: 30242133; array job: 30242143 (`0-377%50`).
- Raw root: `/vast/blangme2/jbonnie/hammock/claude-ref-comparison/results/exp_a_rehash_selector64`.
- Figure-cell CSV SHA-256: `5306cd2490c9d584349963a869528d36b3f4f6c9a5c6eea4ecbd031dbff0a153`.
- Expanded seed-42 table SHA-256: `ba50f2bf98474dea0bf66a825f802fbb247ec758dba52280e3c8e8aafa437752`.
- Final Figure 5 SHA-256: `61e00a0f9dc2499b46bdf27d695aa28be9f04ea7acd50c88f4d07c49cf2079f6`.

`input_sha256.txt` records all 18 FASTA checksums. Curated per-seed and
cross-seed results are in `results/`; full matrices and logs remain in the raw
root above.
