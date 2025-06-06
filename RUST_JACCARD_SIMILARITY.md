# Jaccard Similarity in Rust HyperLogLog (HLL)

This document describes how Jaccard similarity is calculated in the Rust HLL implementation and how to adjust the calculation for your specific use case.

## Background

Jaccard similarity is a measure of similarity between two sets, defined as the size of their intersection divided by the size of their union:

```
J(A, B) = |A ∩ B| / |A ∪ B|
```

When using probabilistic data structures like HyperLogLog, calculating exact Jaccard similarity is not directly possible. Instead, we estimate it using the registers stored in the HLL sketches.

## Accuracy Challenges

HyperLogLog-based Jaccard estimation has some inherent limitations:

1. **Compressed Value Range**: Raw Jaccard estimates from HLL tend to be compressed into a narrow range (often around 0.2-0.4), making it difficult to distinguish between truly high and low similarities.

2. **Set Size Sensitivity**: Accuracy varies based on the relative sizes of the two sets being compared. Sets of very different sizes present particular challenges.

3. **Precision Impact**: Higher precision settings (more registers) generally improve cardinality estimates but can still yield poor Jaccard estimates.

## Available Methods

Our implementation provides several methods for Jaccard similarity estimation:

1. **Standard (`jaccard`)**: The original method using min/max of registers to construct intersection and union sketches.

2. **Register Match (`jaccard_register_match`)**: Compares register values directly to estimate similarity - counts matching register values among active registers.

3. **MinMax (`jaccard_minmax`)**: Combines register matching with cardinality ratio adjustment to handle sets of different sizes better.

4. **Improved (`jaccard_improved`)**: Our recommended method that adapts its approach based on set characteristics, using different strategies based on the size ratio of the sets.

5. **MinHash-inspired (`jaccard_with_method(other, "minhash")`)**: Uses register values as hash values in a MinHash-like approach.

## Accessing Different Methods

The Rust implementation provides several ways to access different Jaccard methods:

### Using the Rust Module Directly

```python
import rust_hll

# Create sketches
sketch1 = rust_hll.RustHLL(precision=14)
sketch2 = rust_hll.RustHLL(precision=14)

# Add data to sketches
# ...

# Standard method (min/max registers)
jaccard_standard = sketch1.jaccard(sketch2)

# Register matching method
jaccard_register = sketch1.jaccard_register_match(sketch2)

# MinMax method (register match with cardinality adjustment)
jaccard_minmax = sketch1.jaccard_minmax(sketch2)

# Improved method (adaptive approach)
jaccard_improved = sketch1.jaccard_improved(sketch2)

# MinHash-inspired method
jaccard_minhash = sketch1.jaccard_with_method(sketch2, "minhash")

# Generic method selector
jaccard_any = sketch1.jaccard_with_method(sketch2, "standard")  # or "register_match", "minmax", "improved", "minhash"
```

### Using the Python Wrapper

The Python wrapper (`RustHLL`) primarily uses the improved method with optional empirical correction:

```python
from hammock.lib.rusthll import RustHLL

# Default behavior with correction
sketch1 = RustHLL(precision=14)
sketch2 = RustHLL(precision=14)

# Add data to sketches
# ...

# Get Jaccard similarity with default correction (uses jaccard_improved internally)
jaccard = sketch1.jaccard(sketch2)
# or equivalently:
jaccard = sketch1.estimate_jaccard(sketch2)

# Disable the empirical correction
sketch1_no_corr = RustHLL(precision=14, apply_jaccard_correction=False)
sketch2_no_corr = RustHLL(precision=14, apply_jaccard_correction=False)

# Add the same data
# ...

# Get raw Jaccard value without Python-side correction
raw_jaccard = sketch1_no_corr.jaccard(sketch2_no_corr)
```

**Note**: The Python wrapper (`RustHLL`) does not directly expose all the individual Jaccard methods (`jaccard_register_match`, `jaccard_minmax`, etc.). To access these methods, you need to use the Rust module directly (`rust_hll.RustHLL`) or access the underlying Rust sketch:

```python
from hammock.lib.rusthll import RustHLL

sketch1 = RustHLL(precision=14)
sketch2 = RustHLL(precision=14)

# Add data...

# Access the underlying Rust sketch for more methods
if sketch1.is_using_rust() and sketch2.is_using_rust():
    rust_sketch1 = sketch1._rust_sketch
    rust_sketch2 = sketch2._rust_sketch
    
    # Now you can use all the Jaccard methods
    register_jaccard = rust_sketch1.jaccard_register_match(rust_sketch2)
    minmax_jaccard = rust_sketch1.jaccard_minmax(rust_sketch2)
```

## Method Descriptions

### Improved Method (`jaccard_improved`)

This is the default method used by the Python wrapper. It adapts its strategy based on the cardinality ratio of the two sets:

- For sets with similar sizes (ratio > 0.8): Uses the traditional HLL Jaccard (min/max registers)
- For moderately different sizes (ratio > 0.5): Uses a weighted combination of 60% HLL Jaccard + 40% register matching
- For very different sizes (ratio ≤ 0.5): Uses a weighted combination of 30% HLL Jaccard + 70% register matching

### Register Match Method (`jaccard_register_match`)

This method counts matching register values among active registers (registers where at least one sketch has a non-zero value). It tends to work better for sets of very different sizes.

### MinMax Method (`jaccard_minmax`)

This method combines the register matching approach with a cardinality ratio adjustment, scaling the register Jaccard by `(0.5 + 0.5 * cardinality_ratio)` to account for size differences.

## Choosing the Right Approach

Consider the following guidelines:

1. **For sets of similar size with high similarity**: The improved method with correction generally works well.

2. **For sets with very different sizes**: The raw register match or minmax methods may be more reliable, but expect compressed values.

3. **For critical applications**: When accuracy is crucial, consider using a different similarity measure or non-probabilistic method if possible.

4. **For experimentation**: Use the Rust module directly to access all methods and compare their performance on your specific data.

## Known Limitations

- Jaccard values for sets with low similarity (< 0.2) tend to be overestimated.
- Jaccard values for sets with high similarity (> 0.7) can be underestimated.
- Sets with very different sizes produce less reliable estimates.
- The Python wrapper's empirical correction is tuned for general use cases and may not be optimal for all scenarios.

## Troubleshooting

If you find that Jaccard values are consistently:

- **Too high**: Use `apply_jaccard_correction=False` to disable the empirical correction, or try the register match method directly.
- **Too low**: Double-check your sketch precision and hash size settings. Higher precision can help but has diminishing returns. Try the standard or improved methods.
- **Inconsistent**: Consider if your data characteristics are changing significantly between comparisons. Try different methods to see which performs best for your use case.

## Further Improvements

We continue to refine our Jaccard similarity implementation. If you have suggestions or encounter issues, please report them to help us improve the accuracy for diverse use cases. 