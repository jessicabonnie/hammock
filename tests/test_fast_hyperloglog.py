#!/usr/bin/env python3
"""
Tests for FastHyperLogLog with C++/Cython acceleration.

This module tests the FastHyperLogLog implementation that provides 2-15x performance 
improvements over the standard HyperLogLog through optional C++ acceleration with
Cython fallback.
"""

import pytest
import numpy as np
import time
from typing import List
import warnings

# Import the classes to test
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.hyperloglog_fast import (
    FastHyperLogLog, 
    create_fast_hyperloglog, 
    get_performance_info,
    CPP_AVAILABLE,
    CYTHON_AVAILABLE
)


@pytest.fixture
def test_data_small():
    """Small test dataset for quick tests."""
    return [f"test_item_{i}" for i in range(100)]


@pytest.fixture  
def test_data_medium():
    """Medium test dataset for performance tests."""
    return [f"test_item_{i}" for i in range(1000)]


@pytest.fixture
def test_data_large():
    """Large test dataset for stress tests."""
    return [f"test_item_{i}" for i in range(10000)]


@pytest.mark.quick
class TestFastHyperLogLogBasic:
    """Basic functionality tests for FastHyperLogLog."""
    
    def test_initialization(self):
        """Test basic initialization of FastHyperLogLog."""
        # Test default initialization
        sketch = FastHyperLogLog()
        assert sketch.precision == 8
        assert sketch.num_registers == 256
        assert sketch.kmer_size == 0
        assert sketch.hash_size == 64  # Updated default for better C++ performance
        
        # Test custom initialization
        sketch = FastHyperLogLog(precision=12, kmer_size=21, hash_size=64)
        assert sketch.precision == 12
        assert sketch.num_registers == 4096
        assert sketch.kmer_size == 21
        assert sketch.hash_size == 64
    
    def test_add_string_basic(self, test_data_small):
        """Test basic add_string functionality."""
        sketch = FastHyperLogLog(precision=12)
        
        # Add strings individually
        for item in test_data_small:
            sketch.add_string(item)
        
        cardinality = sketch.estimate_cardinality()
        
        # Should estimate around 100 with reasonable error
        error = abs(cardinality - len(test_data_small)) / len(test_data_small)
        assert error < 0.2, f"Error too high: {error:.3f}"
    
    def test_add_batch_basic(self, test_data_medium):
        """Test basic add_batch functionality."""
        sketch = FastHyperLogLog(precision=12)
        
        # Add all strings in one batch
        sketch.add_batch(test_data_medium)
        
        cardinality = sketch.estimate_cardinality()
        
        # Should estimate around 1000 with reasonable error
        error = abs(cardinality - len(test_data_medium)) / len(test_data_medium)
        assert error < 0.15, f"Error too high: {error:.3f}"
    
    def test_add_kmers_batch(self):
        """Test add_kmers_batch functionality."""
        sketch = FastHyperLogLog(precision=12, kmer_size=0)  # Use whole string mode for direct k-mer processing
        
        # Test k-mer batch processing
        kmers = ["AAA", "AAT", "ATG", "GCC", "TTT"]
        sketch.add_kmers_batch(kmers)
        
        cardinality = sketch.estimate_cardinality()
        
        # Should estimate around 5 k-mers
        error = abs(cardinality - len(kmers)) / len(kmers)
        assert error < 0.3, f"Error too high: {error:.3f}"
    
    def test_cpp_integration(self):
        """Test that C++ integration works when available."""
        if CPP_AVAILABLE:
            # Test C++ acceleration explicitly
            sketch = FastHyperLogLog(precision=12, use_cpp=True)
            assert sketch._acceleration_type == 'C++'
            
            # Test basic functionality with C++
            test_data = ["item1", "item2", "item3", "item4", "item5"]
            sketch.add_batch(test_data)
            cardinality = sketch.estimate_cardinality()
            
            # Should estimate around 5 items (C++ may use different correction)
            error = abs(cardinality - len(test_data)) / len(test_data)
            assert error < 1.0, f"Error too high: {error:.3f} (C++ cardinality: {cardinality:.2f})"
        else:
            # Skip test if C++ not available
            pytest.skip("C++ extension not available")
    
    def test_acceleration_fallback(self):
        """Test that acceleration fallback works correctly."""
        # Test auto-detection
        sketch = FastHyperLogLog(precision=12)
        
        # Should use best available acceleration
        if CPP_AVAILABLE:
            assert sketch._acceleration_type == 'C++'
        elif CYTHON_AVAILABLE:
            assert sketch._acceleration_type == 'Cython'
        else:
            assert sketch._acceleration_type == 'Python'
    
    def test_empty_sketch(self):
        """Test behavior with empty sketch."""
        sketch = FastHyperLogLog(precision=12)
        
        cardinality = sketch.estimate_cardinality()
        assert cardinality == 0.0, f"Empty sketch should return 0, got {cardinality}"
        
        # Test empty batch
        sketch.add_batch([])
        cardinality = sketch.estimate_cardinality()
        assert cardinality == 0.0, f"Empty batch should keep cardinality at 0, got {cardinality}"
    
    def test_performance_comparison(self, test_data_medium):
        """Test that accelerated version is faster than pure Python."""
        if CPP_AVAILABLE or CYTHON_AVAILABLE:
            # Test accelerated version
            start_time = time.time()
            accel_sketch = FastHyperLogLog(precision=12)
            accel_sketch.add_batch(test_data_medium)
            accel_time = time.time() - start_time
            accel_cardinality = accel_sketch.estimate_cardinality()
            
            # Test pure Python version
            start_time = time.time()
            python_sketch = FastHyperLogLog(precision=12, use_cpp=False, use_cython=False)
            python_sketch.add_batch(test_data_medium)
            python_time = time.time() - start_time
            python_cardinality = python_sketch.estimate_cardinality()
            
            # Accelerated version should be faster
            assert accel_time < python_time, f"Accelerated version ({accel_time:.4f}s) should be faster than Python ({python_time:.4f}s)"
            
            # Results should be similar (within 10% error)
            error = abs(accel_cardinality - python_cardinality) / python_cardinality
            assert error < 0.1, f"Results too different: {error:.3f}"
        else:
            pytest.skip("No acceleration available for comparison")


class TestAccelerationIntegration:
    """Tests for C++/Cython acceleration integration."""
    
    def test_acceleration_availability_detection(self):
        """Test that acceleration availability is correctly detected."""
        info = get_performance_info()
        
        assert 'cpp_available' in info
        assert 'cython_available' in info
        assert 'performance_gain' in info
        assert isinstance(info['cpp_available'], bool)
        assert isinstance(info['cython_available'], bool)
        
        if info['cpp_available']:
            assert '2-15x' in info['performance_gain']
        elif info['cython_available']:
            assert '2-5x' in info['performance_gain']
        else:
            assert 'Install' in info['performance_gain']
    
    def test_explicit_acceleration_control(self, test_data_small):
        """Test explicit control over acceleration usage."""
        # Test auto-detection (default)
        sketch_auto = FastHyperLogLog(precision=12)
        sketch_auto.add_batch(test_data_small)
        cardinality_auto = sketch_auto.estimate_cardinality()
        
        # Test explicit Python mode
        sketch_python = FastHyperLogLog(precision=12, use_cpp=False, use_cython=False)
        sketch_python.add_batch(test_data_small)
        cardinality_python = sketch_python.estimate_cardinality()
        
        # Results should be very similar (same seed, same algorithm)
        relative_diff = abs(cardinality_auto - cardinality_python) / cardinality_python
        assert relative_diff < 0.05, f"Auto vs Python results differ by {relative_diff:.3f}"
        
        # Test explicit C++ mode (if available)
        if CPP_AVAILABLE:
            sketch_cpp = FastHyperLogLog(precision=12, use_cpp=True)
            sketch_cpp.add_batch(test_data_small)
            cardinality_cpp = sketch_cpp.estimate_cardinality()
            
            # C++ and Python should give similar results
            relative_diff = abs(cardinality_cpp - cardinality_python) / cardinality_python
            assert relative_diff < 0.1, f"C++ vs Python results differ by {relative_diff:.3f}"
        else:
            # Should raise ImportError when forcing C++ if not available
            with pytest.raises(ImportError, match="C\\+\\+ acceleration explicitly requested"):
                FastHyperLogLog(precision=12, use_cpp=True)
        
        # Test explicit Cython mode (if available)
        if CYTHON_AVAILABLE:
            sketch_cython = FastHyperLogLog(precision=12, use_cython=True)
            sketch_cython.add_batch(test_data_small)
            cardinality_cython = sketch_cython.estimate_cardinality()
            
            # Cython and Python should give identical results (same seed, implementation)
            relative_diff = abs(cardinality_cython - cardinality_python) / cardinality_python
            assert relative_diff < 0.01, f"Cython vs Python results differ by {relative_diff:.3f}"
        else:
            # Should raise ImportError when forcing Cython if not available
            with pytest.raises(ImportError, match="Cython acceleration explicitly requested"):
                FastHyperLogLog(precision=12, use_cython=True)
    
    @pytest.mark.skipif(not (CPP_AVAILABLE or CYTHON_AVAILABLE), reason="No acceleration available")
    def test_performance_improvement(self, test_data_large):
        """Test that acceleration actually improves performance."""
        # Measure Python implementation time
        sketch_python = FastHyperLogLog(precision=12, use_cpp=False, use_cython=False)
        
        start_time = time.time()
        sketch_python.add_batch(test_data_large)
        python_time = time.time() - start_time
        
        # Measure accelerated implementation time
        sketch_accel = FastHyperLogLog(precision=12)
        
        start_time = time.time()
        sketch_accel.add_batch(test_data_large)
        accel_time = time.time() - start_time
        
        # Accelerated should be faster
        speedup = python_time / accel_time
        print(f"Performance test: Python={python_time:.4f}s, {sketch_accel._acceleration_type}={accel_time:.4f}s, Speedup={speedup:.2f}x")
        
        # For large datasets, we expect at least some speedup
        assert speedup > 1.2, f"Expected speedup > 1.2x, got {speedup:.2f}x"
    
    def test_performance_info_method(self, test_data_small):
        """Test the get_performance_info method on instances."""
        sketch = FastHyperLogLog(precision=12)
        sketch.add_batch(test_data_small)
        
        perf_info = sketch.get_performance_info()
        
        assert 'acceleration_type' in perf_info
        assert 'cpp_available' in perf_info
        assert 'cython_available' in perf_info
        assert 'precision' in perf_info
        assert 'expected_speedup' in perf_info
        
        # Check data consistency
        assert perf_info['precision'] == 12
        assert perf_info['cpp_available'] == CPP_AVAILABLE
        assert perf_info['cython_available'] == CYTHON_AVAILABLE
    
    @pytest.mark.skipif(not (CPP_AVAILABLE or CYTHON_AVAILABLE), reason="No acceleration available")
    def test_built_in_benchmarking(self, test_data_medium):
        """Test the built-in benchmarking functionality."""
        sketch = FastHyperLogLog(precision=12)
        
        results = sketch.benchmark_performance(test_data_medium)
        
        # Check benchmark results structure
        expected_keys = ['accelerated_time', 'python_time', 'speedup', 'accelerated_cardinality', 
                        'python_cardinality', 'accuracy_error_percent', 'num_strings', 'acceleration_used']
        
        for key in expected_keys:
            assert key in results, f"Missing key in benchmark results: {key}"
        
        # Check benchmark results validity
        assert results['speedup'] > 0
        assert results['accuracy_error_percent'] < 5.0  # Less than 5% error
        assert results['num_strings'] == len(test_data_medium)
        assert results['acceleration_used'] in ['C++', 'Cython']


class TestCompatibility:
    """Tests for compatibility with existing HyperLogLog."""
    
    def test_interface_compatibility(self, test_data_small):
        """Test that FastHyperLogLog has compatible interface with HyperLogLog."""
        fast_sketch = FastHyperLogLog(precision=12, seed=42)
        regular_sketch = HyperLogLog(precision=12, seed=42)
        
        # Both should have the same basic methods
        methods = ['add_string', 'estimate_cardinality', 'merge', 'is_empty']
        for method in methods:
            assert hasattr(fast_sketch, method), f"FastHyperLogLog missing method: {method}"
            assert hasattr(regular_sketch, method), f"HyperLogLog missing method: {method}"
    
    def test_cardinality_estimation_consistency(self, test_data_medium):
        """Test that cardinality estimates are consistent between implementations."""
        fast_sketch = FastHyperLogLog(precision=14, seed=42, use_cython=False)
        regular_sketch = HyperLogLog(precision=14, seed=42)
        
        # Add same data to both
        fast_sketch.add_batch(test_data_medium)
        for item in test_data_medium:
            regular_sketch.add_string(item)
        
        fast_cardinality = fast_sketch.estimate_cardinality()
        regular_cardinality = regular_sketch.estimate_cardinality()
        
        # Should be very similar (same precision, seed, algorithm)
        relative_diff = abs(fast_cardinality - regular_cardinality) / regular_cardinality
        assert relative_diff < 0.02, f"FastHyperLogLog and HyperLogLog results differ by {relative_diff:.3f}"
    
    def test_merge_compatibility(self, test_data_small):
        """Test that merging works between FastHyperLogLog and HyperLogLog."""
        fast_sketch = FastHyperLogLog(precision=12, seed=42, use_cython=False)
        regular_sketch = HyperLogLog(precision=12, seed=42)
        
        # Add different data to each
        fast_sketch.add_batch(test_data_small[:50])
        for item in test_data_small[50:]:
            regular_sketch.add_string(item)
        
        # Merge regular into fast
        merged_sketch = fast_sketch.merge_new(regular_sketch)
        merged_cardinality = merged_sketch.estimate_cardinality()
        
        # Should estimate around 100 unique items
        error = abs(merged_cardinality - len(test_data_small)) / len(test_data_small)
        assert error < 0.2, f"Merge error too high: {error:.3f}"
    
    def test_similarity_functions(self, test_data_small):
        """Test similarity functions work correctly."""
        sketch1 = FastHyperLogLog(precision=14, seed=42, use_cython=False)
        sketch2 = FastHyperLogLog(precision=14, seed=42, use_cython=False)
        
        # Identical data - should have high Jaccard similarity
        sketch1.add_batch(test_data_small)
        sketch2.add_batch(test_data_small)
        
        jaccard = sketch1.estimate_jaccard(sketch2)
        
        # Should be close to 1.0 for identical sets
        assert jaccard > 0.8, f"Jaccard for identical sets should be high, got {jaccard:.3f}"
        
        # Test with different data
        sketch3 = FastHyperLogLog(precision=14, seed=42, use_cython=False)
        different_data = [f"different_item_{i}" for i in range(100)]
        sketch3.add_batch(different_data)
        
        jaccard_different = sketch1.estimate_jaccard(sketch3)
        
        # Should be close to 0.0 for disjoint sets
        assert jaccard_different < 0.2, f"Jaccard for disjoint sets should be low, got {jaccard_different:.3f}"


class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    def test_create_fast_hyperloglog(self, test_data_small):
        """Test the create_fast_hyperloglog convenience function."""
        # Test basic creation
        sketch = create_fast_hyperloglog(precision=12)
        assert isinstance(sketch, FastHyperLogLog)
        assert sketch.precision == 12
        
        # Test with additional parameters (use whole string mode)
        sketch = create_fast_hyperloglog(precision=14, kmer_size=0, debug=True)
        assert sketch.precision == 14
        assert sketch.kmer_size == 0
        
        # Test functionality
        sketch.add_batch(test_data_small)
        cardinality = sketch.estimate_cardinality()
        error = abs(cardinality - len(test_data_small)) / len(test_data_small)
        assert error < 0.2, f"Error too high: {error:.3f}"
    
    def test_global_performance_info(self):
        """Test the global get_performance_info function."""
        info = get_performance_info()
        
        required_keys = ['cython_available', 'recommended_class', 'performance_gain']
        for key in required_keys:
            assert key in info, f"Missing key in performance info: {key}"
        
        assert info['recommended_class'] == 'FastHyperLogLog'
        
        if CYTHON_AVAILABLE:
            assert '2-5x speedup' in info['performance_gain']
        else:
            assert 'Install Cython' in info['performance_gain']


class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_invalid_parameters(self):
        """Test handling of invalid parameters."""
        # Invalid precision
        with pytest.raises(ValueError):
            FastHyperLogLog(precision=2)  # Too low
        
        with pytest.raises(ValueError):
            FastHyperLogLog(precision=65, hash_size=64)  # Too high for hash size
        
        # Invalid k-mer size
        with pytest.raises(ValueError):
            FastHyperLogLog(kmer_size=-1)
        
        # Invalid window size
        with pytest.raises(ValueError):
            FastHyperLogLog(kmer_size=21, window_size=10)  # Window < kmer
    
    def test_large_batch_processing(self):
        """Test processing very large batches."""
        sketch = FastHyperLogLog(precision=16)  # Higher precision for accuracy
        
        # Create large dataset
        large_data = [f"large_item_{i}" for i in range(50000)]
        
        # Should handle large batch without issues
        sketch.add_batch(large_data)
        cardinality = sketch.estimate_cardinality()
        
        # Should estimate reasonably close to 50000
        error = abs(cardinality - 50000) / 50000
        assert error < 0.1, f"Large batch error too high: {error:.3f}"
    
    def test_memory_efficiency(self, test_data_large):
        """Test that memory usage is reasonable."""
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Create multiple sketches and process data
        sketches = []
        for i in range(10):
            sketch = FastHyperLogLog(precision=12)
            sketch.add_batch(test_data_large)
            sketches.append(sketch)
        
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_increase = final_memory - initial_memory
        
        # Memory increase should be reasonable (less than 100MB for 10 sketches)
        assert memory_increase < 100, f"Memory increase too high: {memory_increase:.1f}MB"
    
    def test_string_handling_edge_cases(self):
        """Test handling of various string types and edge cases."""
        sketch = FastHyperLogLog(precision=12)
        
        # Test various string types
        edge_case_strings = [
            "",           # Empty string
            "a",          # Single character  
            "🎉🚀",        # Unicode/emoji
            "a" * 1000,   # Very long string
            "test\0null", # String with null byte
            "test\nline", # String with newline
        ]
        
        # Should handle all edge case strings without error
        sketch.add_batch(edge_case_strings)
        cardinality = sketch.estimate_cardinality()
        
        # Should estimate around the number of unique strings
        assert cardinality > 0, "Should estimate > 0 for non-empty batch"


class TestIntegrationWithIntervals:
    """Tests for integration with IntervalSketch."""
    
    def test_intervalsketch_uses_fasthll(self):
        """Test that IntervalSketch automatically uses FastHyperLogLog."""
        from hammock.lib.intervals import IntervalSketch
        
        # Create IntervalSketch - should use FastHyperLogLog by default
        sketch = IntervalSketch(mode='A', precision=12)
        
        # Check that it's using FastHyperLogLog
        assert isinstance(sketch.sketch, FastHyperLogLog), \
            f"Expected FastHyperLogLog, got {type(sketch.sketch)}"
    
    def test_intervalsketch_fallback_control(self):
        """Test controlling FastHyperLogLog usage in IntervalSketch."""
        from hammock.lib.intervals import IntervalSketch
        
        # Test explicit disable
        sketch = IntervalSketch(mode='A', precision=12, use_fast_hll=False)
        assert isinstance(sketch.sketch, HyperLogLog), \
            f"Expected HyperLogLog, got {type(sketch.sketch)}"
        assert not isinstance(sketch.sketch, FastHyperLogLog), \
            "Should not use FastHyperLogLog when explicitly disabled"


@pytest.mark.full
class TestPerformanceBenchmarks:
    """Slow performance benchmark tests."""
    
    @pytest.mark.skipif(not CYTHON_AVAILABLE, reason="Cython not available")
    def test_comprehensive_performance_benchmark(self):
        """Comprehensive performance benchmark comparing implementations."""
        test_sizes = [100, 1000, 5000, 10000]
        results = []
        
        for size in test_sizes:
            test_data = [f"benchmark_item_{i}" for i in range(size)]
            
            # Benchmark Python implementation
            sketch_python = FastHyperLogLog(precision=14, use_cython=False)
            start_time = time.time()
            sketch_python.add_batch(test_data)
            python_time = time.time() - start_time
            
            # Benchmark Cython implementation  
            sketch_cython = FastHyperLogLog(precision=14, use_cython=True)
            start_time = time.time()
            sketch_cython.add_batch(test_data)
            cython_time = time.time() - start_time
            
            speedup = python_time / cython_time if cython_time > 0 else float('inf')
            
            results.append({
                'size': size,
                'python_time': python_time,
                'cython_time': cython_time,
                'speedup': speedup
            })
            
            print(f"Size {size}: Python={python_time:.4f}s, Cython={cython_time:.4f}s, Speedup={speedup:.2f}x")
        
        # Verify performance improves with size
        large_speedup = next(r['speedup'] for r in results if r['size'] == 10000)
        assert large_speedup > 1.5, f"Expected significant speedup for large datasets, got {large_speedup:.2f}x"


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 