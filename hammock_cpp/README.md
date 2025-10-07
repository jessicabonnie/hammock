# Hammock C++ Implementation

A high-performance C++ implementation of the Hammock genomic interval comparison tool using HyperLogLog sketches with OpenMP parallelization. Supports interval-based (Mode A), point-based (Mode B), and combined (Mode C) comparison.

## Quick Start

### Build

```bash
cd /home/jbonnie1/interval_sketch/hammock/hammock_cpp
make
```

### Run Tests

```bash
make test
```

This runs 36 unit tests covering:
- Shared utilities (parsing, normalization)
- Mode A (interval-based comparison)
- Mode B (point-based comparison)
- Mode B subsampling (--subB)
- Mode C (combined comparison)
- Mode C interval expansion (--expA)
- Cross-mode behavior validation

### Run Examples

```bash
# Mode A: Compare intervals
make example_a

# Mode B: Compare points (bases)
make example_b

# Mode C: Compare combined intervals+points
make example_c

# Run all three
make example
```

## Modes

### Mode A: Interval-Based Comparison

Compares genomic intervals as discrete units. Each interval (chr, start, end) is treated as a single element regardless of its size.

**Use when**: You want to know how many intervals are shared between files, independent of their size.

**Example**:
```
chr1:1000-2000  →  "1-1000-2000-A"  (1 element)
chr1:1000-10000 →  "1-1000-10000-A" (1 element)
```

### Mode B: Point-Based Comparison

Compares all individual base positions within intervals. Each base position becomes a separate element.

**Use when**: You want to know what fraction of the genome is covered by both files.

**Example**:
```
chr1:1000-2000  →  1000, 1001, 1002, ..., 1999  (1000 elements)
chr1:1000-10000 →  1000, 1001, 1002, ..., 9999  (9000 elements)
```

**Subsampling (--subB)**: For large intervals, Mode B can be slow. Use `--subB` to randomly sample points:
- `--subB 1.0`: Process all points (default)
- `--subB 0.5`: Process 50% of points (faster)
- `--subB 0.1`: Process 10% of points (much faster)

**Key property**: Subsampling is **deterministic** - the same point is always selected or rejected across different files, preserving meaningful Jaccard similarities.

### Mode C: Combined Comparison

Compares both intervals AND points in a single sketch. Each interval becomes one element, and each base position becomes another element.

**Use when**: You want a similarity measure that considers both structural overlap (intervals) and base-level coverage (points).

**Example**:
```
chr1:1000-2000  →  "1-1000-2000-A" (1 interval)
                   + 1000, 1001, ..., 1999 (1000 points)
                   = 1001 total elements
```

**Behavior**: Mode C combines the best of both worlds - it captures interval-level similarity while also accounting for the actual genomic coverage.

**Interval Expansion (--expA)**: By default, each interval contributes 1 element, which can be dominated by points when intervals are large. Use `--expA` to expand each interval:
- `--expA 0`: 1 version per interval (default)
- `--expA 0.5`: 3.16 versions per interval (10^0.5)
- `--expA 1.0`: 10 versions per interval (10^1)
- `--expA 1.5`: 31.6 versions per interval (10^1.5)
- `--expA 2.0`: 100 versions per interval (10^2)

**Note**: `--expA` accepts decimal values for fine-grained control over interval weighting.

This allows you to balance the relative weight of intervals vs. points in the combined sketch.

### Comparison

Using the example data (comparing example_test2.bed vs example_primary.bed):

| Mode | Cardinality | Jaccard | Meaning | Speed |
|------|-------------|---------|---------|-------|
| A    | ~8, ~7      | 0.36    | 36% of intervals are shared | Fast |
| B (subB=1.0) | ~8K, ~7K | 0.41 | 41% of bases are shared | Slower |
| B (subB=0.1) | ~800, ~700 | ~0.38 | ~38% of bases are shared | Much faster |
| C (default) | ~8K, ~7K | 0.41 | Combined similarity (point-dominated) | Slower |
| C (expA=2) | ~8.8K, ~7.7K | 0.42 | Balanced interval+point weight | Slower |
| C (expA=2, subB=0.1) | ~1.6K, ~1.4K | 0.37 | Balanced + fast | Much faster |

**Notes**: 
- Mode C cardinality is dominated by points unless `--expA` is used
- `--expA 2` adds 100 versions of each interval to balance with point contributions
- Mode C Jaccard reflects both interval and point overlap
- Subsampling slightly changes Jaccard but preserves relative similarities

## Usage

```bash
hammock <files_list> <primary_list> [options]

Options:
  --mode <A|B|C>         : Comparison mode (default: A)
                           A = interval-based (chr-start-end)
                           B = point-based (all bases in intervals)
                           C = combined (intervals + points)
  --subB <float>         : Subsampling rate for points (0.0-1.0, default: 1.0)
                           Used in Mode B and Mode C
  --expA <float>         : Interval expansion exponent (default: 0.0)
                           Adds 10^expA versions of each interval (Mode C only)
                           Accepts decimal values (e.g., 0.5, 1.5, 2.5)
                           Balances interval weight relative to points
  --precision, -p <int>  : HyperLogLog precision (default: 18)
  --separator, -s <str>  : Separator for strings (default: "-")
  --output, -o <prefix>  : Output file prefix (default: hammock)
                           Naming: {prefix}_hll_p{precision}_jacc{mode}[_expA{expA}][_B{subB}].csv
                           Matches Python hammock naming convention
  --verbose, -v          : Print verbose progress information
  --help, -h             : Show help message

Examples:
  # Mode A: Compare intervals
  ./bin/hammock files.txt primary.txt --mode A -o results_A.csv

  # Mode B: Compare all points
  ./bin/hammock files.txt primary.txt --mode B -o results_B.csv
  
  # Mode B with subsampling (faster)
  ./bin/hammock files.txt primary.txt --mode B --subB 0.1 -o results_B.csv
  
  # Mode C: Combined intervals and points
  ./bin/hammock files.txt primary.txt --mode C -o results_C.csv
  
  # Mode C with interval expansion (balanced weighting)
  ./bin/hammock files.txt primary.txt --mode C --expA 2.0 -o results
    → results_hll_p18_jaccC_expA2.00.csv
  
  # Mode C with decimal expansion for fine-tuning
  ./bin/hammock files.txt primary.txt --mode C --expA 1.5 -o results
    → results_hll_p18_jaccC_expA1.50.csv
  
  # Mode C with expansion and subsampling (balanced and fast)
  ./bin/hammock files.txt primary.txt --mode C --expA 1.5 --subB 0.25 -o results
    → results_hll_p18_jaccC_expA1.50_B0.25.csv
```

## Code Architecture

### Unified Design

The codebase uses a **unified architecture** with shared utilities and mode-specific processing:

```
hammock.cpp
├── Shared Utilities (used by both modes)
│   ├── normalize_chromosome()      # Remove "chr" prefix
│   ├── is_header_or_blank()        # Filter headers
│   ├── parse_bed_line()            # Parse BED format
│   ├── read_filepath_list()        # Read file lists
│   └── calculate_jaccard()         # Compute similarity
│
├── Mode A: Interval-Based
│   ├── create_interval_string()    # Format: chr-start-end-A
│   └── process_bed_file_mode_a()   # Hash intervals
│
├── Mode B: Point-Based
│   ├── create_point_string()       # Format: chr-pos-B
│   └── process_bed_file_mode_b()   # Hash each base position
│
├── Mode C: Combined
│   └── process_bed_file_mode_c()   # Hash intervals + points
│
└── Main Program
    └── Unified CLI handling all three modes
```

### Why Unified?

**Code reuse benefits**:
- ~200 lines of shared utility code (60% of codebase)
- Single binary supports both modes
- Consistent behavior across modes
- Easier maintenance and testing

**Avoided duplication**:
- BED file parsing
- Chromosome normalization
- File list reading  
- Jaccard calculation
- Command-line interface
- CSV output generation

**Mode-specific code** (~100 lines each):
- String formatting (interval vs. point)
- Processing loop (1 hash per interval vs. N hashes per interval)

## Build Targets

```bash
make              # Build hammock binary
make test         # Run 36 unit tests (all passing)
make example_a    # Run Mode A example
make example_b    # Run Mode B example
make example      # Run both examples
make install      # Copy binary to hammock_cpp root
make clean        # Clean build artifacts
make clean-all    # Clean everything including HLL
```

## Output Format

The C++ implementation uses the same naming convention as Python hammock for seamless integration:

```
{prefix}_hll_p{precision}_jacc{mode}[_expA{expA:.2f}][_B{subB:.2f}].csv
```

Examples:
- `results_hll_p14_jaccA.csv` (Mode A)
- `results_hll_p14_jaccB_B0.10.csv` (Mode B with 10% subsampling)
- `results_hll_p14_jaccC_expA1.50_B0.25.csv` (Mode C with expA=1.5 and subB=0.25)

CSV columns:
- **All modes**: bed1, bed2, sketch_type, mode, precision, jaccard, intersection, union, cardinality1, cardinality2
- **Mode C only**: Additional columns: expA, subB

## Project Structure

```
hammock_cpp/
├── src/
│   └── hammock.cpp          # Unified implementation (Modes A & B)
├── test/
│   └── test_hammock.cpp     # 21 unit tests
├── examples/
│   ├── data/               # Sample BED files
│   ├── files_list.txt      # Input file list
│   ├── primary_list.txt    # Reference file list
│   ├── example_output_A.csv # Mode A results
│   └── example_output_B.csv # Mode B results
├── bin/
│   ├── hammock             # Main binary
│   └── test_hammock        # Test binary
├── docs/                    # Documentation
├── Makefile                # Build system
└── README.md               # This file
```

## Output Format

CSV with columns:

| Column | Description |
|--------|-------------|
| bed1 | Path to first BED file |
| bed2 | Path to second BED file |
| sketch_type | "hyperloglog" |
| mode | "A" (interval) or "B" (point) |
| precision | HyperLogLog precision used |
| jaccard | Jaccard similarity (0-1) |
| intersection | Estimated intersection size |
| union | Estimated union size |
| cardinality1 | Estimated cardinality of file 1 |
| cardinality2 | Estimated cardinality of file 2 |

### Understanding the Output Values

**Important Note on Accuracy:**

The output values use different HLL estimation methods with varying accuracy:

| Metric | Method | Typical Error | Notes |
|--------|--------|---------------|-------|
| **jaccard** | Register-based comparison | **~2%** | **Most accurate - use this!** |
| cardinality1, cardinality2 | Standard HLL | ~1.6% | Individual set sizes |
| union | MAX of registers | ~1.6% | Standard HLL union |
| intersection | MIN of registers | ~3-5% | Less accurate than Jaccard |

**Key Insight:** Due to different estimation methods, `jaccard ≠ intersection/union` in the output. 

**Recommendation:** The `jaccard` value is the most reliable metric. If you need intersection/union sizes, consider deriving them as:
- `intersection ≈ jaccard × union`
- This ensures mathematical consistency and leverages the accurate Jaccard estimate

**Why this approach?**
- Register-based Jaccard directly compares HLL registers (very accurate: R²=0.9999 vs bedtools)
- Traditional intersection via MIN of registers has higher error
- Inclusion-exclusion (`|A| + |B| - |A∪B|`) compounds errors even more (avoided in this implementation)

See `benchmarks/accuracy_comparison.py` to validate accuracy on your data.

## Performance

### Parallelization Features

The C++ implementation uses **OpenMP** for multi-level parallelization:
- **File-level parallelization**: Multiple files processed simultaneously
- **Interval-level parallelization**: Intervals within each file processed in parallel
- **Jaccard calculations**: All pairwise comparisons computed in parallel
- **Automatic scaling**: Uses all available CPU cores by default

**Intelligent Thread Selection:**

The program automatically selects an optimal number of threads based on workload size:
- **Small workloads** (≤4 comparisons): 4 threads
- **Medium workloads** (≤16 comparisons): 8 threads
- **Large workloads** (≤64 comparisons): 16 threads
- **Very large workloads** (>64 comparisons): 24 threads (capped for efficiency)

This prevents thread overhead on small tasks while maximizing performance on large workloads.

**Manual thread control:**
```bash
# Automatic selection (recommended)
./hammock files.txt primary.txt --mode B -p 14 -o results

# Override with specific thread count
OMP_NUM_THREADS=16 ./hammock files.txt primary.txt --mode B -p 14 -o results

# Force sequential execution (debugging)
OMP_NUM_THREADS=1 ./hammock files.txt primary.txt --mode B -p 14 -o results
```

### Benchmark Results (10,000 intervals/file, 48 cores)

**Mode A vs bedtools:**

| Files | bedtools | hammock_cpp | Speedup |
|-------|----------|-------------|---------|
| 2 | 0.50s | 0.07s | **7x faster** |
| 8 | 2.77s | 0.13s | **21x faster** |
| 16 | 7.97s | 0.13s | **61x faster** |

**Mode B/C Parallelization Scaling:**

| Files | Time | Efficiency |
|-------|------|------------|
| 2 | 11.78s | baseline |
| 4 | 11.72s | 2x files, same time ✓ |
| 8 | 11.83s | 4x files, same time ✓ |
| 16 | 12.07s | 8x files, ~same time ✓ |

### Mode A (Interval-based)
- **Speed**: ~1M intervals/second per core
- **Parallelization**: Near-linear scaling with file count
- **Memory**: 2^precision bytes per sketch
- **Scalability**: Handles millions of intervals efficiently

### Mode B (Point-based)
- **Speed**: ~1M points/second (full sampling)
- **Speed with subsampling**: Processing time reduces proportionally
  - subB=0.5: ~2x faster
  - subB=0.1: ~10x faster
- **Parallelization**: Multi-level (file + interval level)
- **Memory**: Same as Mode A (2^precision bytes)
- **Scalability**: Use subsampling for large genomes
- **Note**: Processing time proportional to `total_interval_length × subB`

**Subsampling recommendations**:
- Large genomes (>1M bases): Use subB=0.1 to 0.3
- Medium genomes: Use subB=0.3 to 0.7
- Small datasets: Use subB=1.0 (no subsampling)

### Memory Usage by Precision

| Precision | Memory | Error  | Use Case |
|-----------|--------|--------|----------|
| 12        | 4 KB   | ~6.4%  | Quick tests |
| 14        | 16 KB  | ~3.2%  | Development |
| 16        | 64 KB  | ~1.6%  | Production |
| 18        | 256 KB | ~0.8%  | High accuracy (default) |
| 20        | 1 MB   | ~0.4%  | Maximum accuracy |

## Testing

The test suite includes:

**Shared Function Tests** (7 tests)
- Chromosome normalization
- BED line parsing
- File list reading
- Jaccard calculation

**Mode A Tests** (4 tests)
- Interval string creation
- File processing
- Jaccard comparisons

**Mode B Tests** (7 tests)
- Point string creation
- Point expansion
- Overlapping intervals
- Jaccard comparisons with points

**Mode B Subsampling Tests** (5 tests)
- Reduces point count proportionally
- Deterministic selection (same point always selected/rejected)
- Different subsampling rates produce expected counts
- Consistent results across runs
- Edge cases (subB=0.0 and subB=1.0)

**Mode C Tests** (6 tests)
- Combines intervals and points
- Comparison with Mode A and B
- Jaccard calculations with combined sketches
- Partial overlap scenarios
- Subsampling in Mode C
- Deterministic subsampling

**Mode C expA Tests** (4 tests)
- Interval expansion (expA=0, 1, 2)
- Weight balancing between intervals and points
- Jaccard with expansion
- Combined expansion and subsampling

**Mode Comparison Tests** (3 tests)
- Verify modes produce different results
- Test size sensitivity differences
- Validate mode-specific behavior

All tests use real file I/O and validate against expected HLL estimates.

## Implementation Notes

### Mode A: Intervals as Elements

Each interval is hashed once:
```cpp
"chr1-1000-2000-A" → hash → add to sketch
```

- Fast: O(n) where n = number of intervals
- Size-independent: Small and large intervals treated equally
- Measures "interval overlap"

### Mode B: Points as Elements

Each base position in each interval is hashed:
```cpp
for (pos = start; pos < end; pos++) {
    // Deterministic subsampling
    if (subsample < 1.0) {
        uint32_t point_hash = hash_32("chr1-pos-B");
        if (point_hash > threshold) continue;
    }
    "chr1-pos-B" → hash → add to sketch
}
```

- Slower: O(sum of interval lengths)
- Size-dependent: Large intervals contribute more elements
- Measures "base-pair overlap"
- **Subsampling**: Hash each point to decide inclusion (deterministic across files)

### Mode C: Combined Elements

Each interval is hashed multiple times (based on `expA`), and each point is hashed separately:
```cpp
// Add interval with expansion (Mode A component)
size_t expansions = static_cast<size_t>(pow(10.0, expA));
for (exp = 0; exp < expansions; exp++) {
    "chr1-1000-2000-A-{exp}" → hash → add to sketch
}

// Add points (Mode B component)
for (pos = 1000; pos < 2000; pos++) {
    if (subsample) { /* deterministic subsampling */ }
    "chr1-pos-B" → hash → add to sketch
}
```

- Moderate speed: O(10^expA × intervals + subB × sum of interval lengths)
- Combines both types of information
- Measures "combined structural and coverage similarity"
- **Interval expansion**: Accepts decimal values (e.g., `--expA 1.5` = 31.6x)
- **Use case**: When both interval boundaries and base-level overlap matter, with fine-grained weighting control

### Design Decisions

1. **Unified vs Separate Binaries**: Chose unified for code reuse and maintainability
2. **String Formatting**: Simple string concatenation (can upgrade to xxHash later)
3. **Mode Suffix**: Added "-A" and "-B" suffixes to prevent cross-mode collisions
4. **Processing Verbosity**: Different progress messages for intervals vs. points
5. **Error Handling**: Shared validation for both modes

## Future Enhancements

- Mode D: Sequence comparison with minimizers
- xxHash integration for exact Python compatibility
- Parallel file processing (multi-threaded)
- Compressed file support (.gz, .bgz)
- BigBed format support

## Comparison with Python Hammock

### Similarities
- ✓ Same algorithms (Mode A and B)
- ✓ Compatible output format
- ✓ Same chromosome normalization

### Advantages
- ⚡ 5-10x faster
- 📦 Standalone (no Python dependencies)
- 💾 Lower memory overhead
- 🚀 SIMD-optimized HLL operations

### Current Limitations
- No Mode D support yet
- Simplified hash function (not xxHash)
- No subA (interval subsampling) yet

## Dependencies

- C++17 compiler (g++ or clang++)
- HLL library (included in `../hll/`)
- SIMD support (SSE2 minimum, AVX2/AVX-512 recommended)

## Related

- Python Hammock: `/home/jbonnie1/interval_sketch/hammock/`
- HLL Library: `/home/jbonnie1/interval_sketch/hammock/hll/`
- Documentation: `docs/`

---

**Version**: 0.7.4 (Determinism and accuracy fixes)  
**Tests**: 39/39 passing ✓ (includes critical bug prevention tests)  
**Status**: Production-ready

### Recent Additions (v0.7.4) ⚠️ CRITICAL ACCURACY FIX
- 🐛 **Fixed Non-Determinism in Mode B/C**:
  - **Issue**: Dynamic OpenMP scheduling caused non-deterministic results (same file produced different cardinalities)
  - **Impact**: Accuracy was severely degraded (R²=0.42 instead of ~0.999)
  - **Fix**: Changed from `schedule(dynamic)` to `schedule(static)` for deterministic chunk assignment
  - **Result**: Perfectly deterministic - same input always produces same output
- 🐛 **Fixed HLL Union Calculation**:
  - **Issue**: HLL library's `operator+` uses broken copy constructor, causing segfaults
  - **Fix**: Avoid copying by creating empty sketch and using `+=` twice
  - **Result**: Union calculations now work correctly
  - **Upgrade urgency**: CRITICAL - all Mode B/C results were affected
- ✅ **Added Critical Bug Prevention Tests**:
  - **Determinism test**: Runs same file 3 times, verifies identical cardinality (catches non-deterministic scheduling)
  - **Union/intersection sanity test**: Verifies union >= max(card1, card2), intersection <= min(card1, card2) (catches broken HLL operations)
  - **Identity test**: Compares file to itself, verifies Jaccard=1.0 and card1=card2 (catches race conditions)
  - **Why this matters**: Previous bugs passed all 36 tests - these new tests specifically target the types of bugs we encountered

### Previous Additions (v0.7.3) 🚀 PERFORMANCE FIX
- ⚡ **Improved Parallelization Strategy for Mode B/C**:
  - **Issue**: At 16 files, premature switch to file-level parallelism caused performance regression
  - **Impact**: 16 files ran slower than expected (~1.5 threads per file instead of 24)
  - **Fix**: Adaptive strategy now considers threads-per-file ratio
    - 8-24 files: Only parallelize if we have 3+ threads per file available
    - Otherwise: Process files sequentially to maximize internal parallelism
    - >24 files: Always parallelize (sequential would be too slow)
  - **Result**: 16 files now stay sequential with full 24-thread parallelism per file
  - **Performance**: Significant improvement at 16-24 files, no regression elsewhere

### Previous Additions (v0.7.2) ⚠️ CRITICAL BUG FIX
- 🐛 **Fixed Race Condition in Mode B and Mode C**: 
  - **Issue**: Multiple threads were modifying the same HLL sketch simultaneously without synchronization
  - **Impact**: Corrupted sketches, nonsensical results (union=cardinality, impossible Jaccard values)
  - **Fix**: Thread-local sketches with merge at the end
  - **Result**: Mode B and Mode C now produce correct, accurate results
  - **Performance**: Minimal overhead (~1μs per thread), maintains full parallelism
  - **Upgrade urgency**: HIGH - all Mode B/C users should upgrade immediately

### Previous Additions (v0.7.1)
- ✅ **Intelligent Thread Selection**: Automatically chooses optimal thread count based on workload
  - Small workloads (≤4 comparisons): 4 threads
  - Medium workloads (≤16 comparisons): 8 threads
  - Large workloads (≤64 comparisons): 16 threads
  - Very large workloads (>64): 24 threads (capped for efficiency)
  - Prevents thread overhead on small tasks while maximizing performance on large workloads
  - Respects user override via `OMP_NUM_THREADS`

### Previous Additions (v0.7.0)
- ✅ **OpenMP Parallelization**: Multi-level parallelization for dramatic speedups
  - File-level: Process multiple files simultaneously
  - Interval-level: Parallelize processing within each file
  - Comparison-level: All pairwise Jaccard calculations in parallel
  - **Result**: 16 files process in ~same time as 2 files (Mode B/C)
  - **Result**: Up to 61x faster than bedtools (Mode A, 16 files)

### Previous Additions (v0.6.0)
- ✅ Decimal values for `--expA` (e.g., 1.5 → 31.6x expansion)
- ✅ Python hammock output naming convention
- ✅ Additional CSV columns (expA, subB) for Mode C
- ✅ Full compatibility with Python hammock downstream analysis

### Previous Additions (v0.5.0)
- ✅ `--expA` parameter for interval expansion in Mode C
- ✅ Controllable weighting of intervals vs. points
- ✅ 4 new tests for expA functionality

### Previous Additions (v0.4.0)
- ✅ Mode C implementation (combined intervals + points)
- ✅ Subsampling support in Mode C
- ✅ 6 new tests for Mode C functionality

### Previous Additions (v0.3.0)
- ✅ Mode B subsampling (--subB parameter)
- ✅ Deterministic point selection via FNV hash
- ✅ Speed improvements for large genomes
