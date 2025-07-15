# cProfile Integration for Hammock Benchmarks

This document explains how to use the new cProfile integration with your existing benchmarks.

## Overview

The profiling integration provides three main ways to profile your code:

1. **ProfiledBenchmark**: Context manager for detailed profiling
2. **@profile_function**: Decorator for automatic profiling
3. **BenchmarkProfiler**: Advanced profiler for multiple configurations

## Files

- `profiling_utils.py`: Core profiling utilities
- `benchmark_sketches_profiled.py`: Enhanced sketch benchmarks with profiling
- `benchmark_hll_profiled.py`: Enhanced HyperLogLog benchmarks with profiling
- `profiling_example.py`: Simple example showing all profiling methods

## Quick Start

### 1. Basic Usage with Context Manager

```python
from profiling_utils import ProfiledBenchmark

with ProfiledBenchmark("my_test", output_dir="results") as profiler:
    # Your benchmark code here
    sketch = HyperLogLog(precision=12)
    for i in range(10000):
        sketch.add_string(f"item_{i}")
    result = sketch.estimate_cardinality()

print(f"Wall time: {profiler.wall_time:.4f}s")
print(f"Memory delta: {profiler.memory_delta:.2f} MiB")
```

### 2. Using the Decorator

```python
from profiling_utils import profile_function

@profile_function(save_profile=True, output_dir="results")
def my_benchmark():
    sketch = HyperLogLog(precision=12)
    for i in range(10000):
        sketch.add_string(f"item_{i}")
    return sketch.estimate_cardinality()

result = my_benchmark()
# Profile info automatically added to result if it's a dictionary
```

### 3. Advanced Multi-Configuration Profiling

```python
from profiling_utils import BenchmarkProfiler

profiler = BenchmarkProfiler(output_dir="results")

configurations = [
    {'name': 'precision_12', 'precision': 12, 'num_items': 10000},
    {'name': 'precision_16', 'precision': 16, 'num_items': 10000},
]

def run_test(precision, num_items):
    sketch = HyperLogLog(precision=precision)
    for i in range(num_items):
        sketch.add_string(f"item_{i}")
    return sketch.estimate_cardinality()

results = profiler.run_profiled_benchmark(
    benchmark_func=run_test,
    name="precision_comparison",
    configurations=configurations
)

# Generate reports and plots
profiler.save_summary_report()
profiler.plot_performance_comparison('wall_time')
```

## Running the Examples

### 1. Simple Example

```bash
cd benchmarks
python profiling_example.py
```

This will:
- Run basic profiling examples
- Generate profile files in `results/`
- Create performance comparison plots
- Show summary statistics

### 2. Enhanced Sketch Benchmarks

```bash
cd benchmarks
python benchmark_sketches_profiled.py
```

This will:
- Profile different sketch types (HyperLogLog, MinHash, etc.)
- Test different parameter values
- Generate detailed performance analysis
- Create comparison reports

### 3. HyperLogLog Deep Analysis

```bash
cd benchmarks
python benchmark_hll_profiled.py
```

This will:
- Profile HyperLogLog implementations in detail
- Compare Python vs Rust implementations (if available)
- Analyze performance bottlenecks
- Generate comprehensive reports

## Understanding the Output

### Profile Files

The profiling generates several types of files:

1. **`.prof` files**: Binary profile data for detailed analysis
2. **`_summary.txt` files**: Human-readable profile summaries
3. **`benchmark_summary.csv`**: Tabular summary of all benchmarks
4. **Performance plots**: Visual comparisons of different configurations

### Key Metrics

- **Wall Time**: Total execution time
- **Memory Delta**: Memory usage change during execution
- **Function Stats**: Per-function timing and call counts
- **Top Functions**: Functions consuming the most time

### Profile Analysis

You can analyze saved profile files:

```python
from profiling_utils import analyze_profile_file

analyze_profile_file('results/my_benchmark.prof')
```

This will show:
- Top functions by cumulative time
- Top functions by self time
- Call statistics
- Detailed timing information

## Integration with Existing Benchmarks

### Minimal Integration

To add profiling to existing benchmarks with minimal changes:

```python
# Original code
def my_benchmark():
    # ... benchmark code ...
    return results

# Add profiling with just two lines
from profiling_utils import ProfiledBenchmark

with ProfiledBenchmark("my_benchmark") as profiler:
    results = my_benchmark()
```

### Full Integration

For complete integration, see `benchmark_sketches_profiled.py` which shows:
- How to wrap existing functions
- How to add profiling to test suites
- How to generate comprehensive reports
- How to maintain backward compatibility

## Performance Impact

The profiling overhead is typically:
- **5-10% slowdown** for most operations
- **Minimal memory overhead** for profile data storage
- **No impact** on accuracy of results

## Best Practices

1. **Profile Representative Workloads**: Use realistic data sizes and patterns
2. **Run Multiple Iterations**: Average results across multiple runs
3. **Focus on Hotspots**: Pay attention to functions with high cumulative time
4. **Compare Implementations**: Use profiling to compare different approaches
5. **Save Profile Data**: Keep profile files for later analysis

## Troubleshooting

### Common Issues

1. **Import Errors**: Make sure all dependencies are installed
2. **Permission Errors**: Check write permissions for results directory
3. **Memory Issues**: Use smaller datasets for initial profiling
4. **Plot Display**: Install matplotlib for visualization features

### Getting Help

If you encounter issues:
1. Check the example files for working code
2. Verify your benchmark code works without profiling first
3. Start with simple profiling before using advanced features
4. Check the results directory for error logs

## Advanced Features

### Custom Analysis

You can extend the profiling utilities:

```python
from profiling_utils import ProfiledBenchmark

with ProfiledBenchmark("test") as profiler:
    # Your code here
    pass

# Custom analysis
stats = profiler.get_function_stats('my_function')
for func, data in stats.items():
    print(f"{func}: {data['call_count']} calls, {data['total_time']:.4f}s")
```

### Comparing Profiles

```python
from profiling_utils import compare_profiles

compare_profiles([
    'results/implementation_a.prof',
    'results/implementation_b.prof'
])
```

This will show side-by-side comparisons of different implementations.

## Next Steps

1. Try the examples to get familiar with the profiling tools
2. Integrate profiling into your existing benchmarks
3. Use the insights to optimize your code
4. Share profile reports with your team for collaborative optimization

Happy profiling! 