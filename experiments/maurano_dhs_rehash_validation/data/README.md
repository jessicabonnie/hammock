# Input policy

Biological inputs are not copied here. `config.yaml` points read-only to the
existing Maurano FASTA list, filename key, and BEDTools reference. Run
`scripts/inventory_inputs.py` to write resolved paths, sizes, and SHA-256
checksums to `results/manifests/inputs.csv`.
