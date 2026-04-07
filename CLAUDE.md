# hammock

A research repository for minimizer-based sequence sketching tools and validation experiments.

## Repo structure

```
hammock/
├── src/           # Core sketching library
├── tests/         # Unit and integration tests
└── experiments/   # Independent validation experiments (each has its own CLAUDE.md)
    ├── claude-ref-comparison/
    └── ...        # Other experiments — do not load their context unless working in them
```

## Key commands

```bash
# Run tests
# TODO: fill in your actual test command

# Build
# TODO: fill in your actual build command
```

## Conventions

- All experiment documents are Markdown-first; LaTeX generated via Pandoc later
- Experiments are self-contained — dependencies, configs, and workflows live inside each experiment directory
- Reproducible workflows use Snakemake + Singularity + SLURM

## Finding things

- For experiment-specific data paths, cluster config, and workflow details: read the experiment's own `CLAUDE.md`
- For sketching tool API and parameters: `src/` and its docstrings
- For bibliography: each experiment's `docs/references.bib`
