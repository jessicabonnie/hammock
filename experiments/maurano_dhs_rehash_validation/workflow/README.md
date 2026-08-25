# Workflow boundary

The checked-in Python drivers are the workflow implementation. Generate the
two dry-run tables before constructing any site-specific scheduler wrapper:

```bash
python scripts/extract_exact_features.py --config config.yaml --dry-run
python scripts/run_rehash_sweep.py --config config.yaml \
  --hammock /isolated/env/bin/hammock --source-commit "$(git rev-parse HEAD)" --dry-run
```

The runner intentionally executes only one named job per invocation. This
makes scheduler concurrency, memory, and authorization visible rather than
embedding cluster-specific submission behavior in the experiment. Do not add
or submit an array until the two-sample and one-cell gates pass and the user
authorizes shared computation.
