# rust_hll

Rust implementation of the HyperLogLog algorithm for the hammock package.

## Building the Extension

### Prerequisites

- Rust toolchain (install via [rustup](https://rustup.rs/))
- Python 3.6+
- maturin (`pip install maturin`)

### Development Build

For local development:

```bash
python -m maturin develop
```

This will build and install the extension in development mode, making it available to Python.

### Release Build

To build a wheel for distribution:

```bash
python -m maturin build --release
```

The resulting wheel file will be in the `target/wheels` directory.

## Troubleshooting

### Installation Issues

If you encounter issues during installation, try these steps:

1. Make sure Rust is properly installed:
   ```bash
   rustc --version
   ```

2. Update Rust to the latest version:
   ```bash
   rustup update
   ```

3. Clear any cached or partially built packages:
   ```bash
   rm -rf target/
   ```

4. Make sure maturin is up-to-date:
   ```bash
   pip install --upgrade maturin
   ```

5. Try building with verbose output:
   ```bash
   python -m maturin develop -v
   ```

### Common Issues

- **Missing compiler**: Make sure rustc and cargo are in your PATH
- **Missing Python headers**: Install python-dev or python-devel package for your OS
- **Incompatible versions**: Ensure your Python and Rust versions are compatible with maturin

## Implementation Notes

This extension implements several HyperLogLog cardinality estimation methods:

- Original HyperLogLog algorithm (Flajolet et al.)
- Improved HyperLogLog with bias correction (Ertl)
- Maximum Likelihood Estimation method (Ertl)

The default method is MLE, which provides the best accuracy in most cases.

## Testing

Once built, you can test the extension with:

```bash
python -c "from hammock import FastHyperLogLog; hll = FastHyperLogLog(precision=8); hll.add('test'); print(hll.estimate())"
```

If you see a reasonable estimate (likely 1.0), the extension is working correctly. 