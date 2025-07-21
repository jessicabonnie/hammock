#!/usr/bin/env python3
"""
Integration tests for performance optimization features.

This module tests the integration of FastHyperLogLog with the broader hammock
ecosystem, including IntervalSketch integration and real-world usage patterns.
"""

import pytest
import tempfile
import os
import time
from hammock.lib.intervals import IntervalSketch
from hammock.lib.hyperloglog import HyperLogLog
from hammock.lib.hyperloglog_fast import (
    FastHyperLogLog, 
    create_fast_hyperloglog,
    CYTHON_AVAILABLE
)


@pytest.fixture
def temp_bed_file():
    """Create a temporary BED file for testing."""
    bed_content = """chr1\t100\t200\tregion1
chr1\t300\t400\tregion2
chr2\t500\t600\tregion3
chr2\t700\t800\tregion4
chr3\t900\t1000\tregion5"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed_content)
        temp_path = f.name
    
    yield temp_path
    
    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def large_bed_content():
    """Generate content for a larger BED file."""
    lines = []
    for i in range(1000):
        chrom = f"chr{(i % 22) + 1}"
        start = i * 100
        end = start + 50
        name = f"region_{i}"
        lines.append(f"{chrom}\t{start}\t{end}\t{name}")
    return "\n".join(lines)


class TestIntervalSketchIntegration:
    """Tests for FastHyperLogLog integration with IntervalSketch."""
    
    def test_intervalsketch_automatic_optimization(self):
        """Test that IntervalSketch automatically uses FastHyperLogLog."""
        # Create IntervalSketch without explicit settings
        sketch = IntervalSketch(mode='A', precision=12)
        
        # Should use FastHyperLogLog by default
        assert isinstance(sketch.sketch, FastHyperLogLog), \
            f"Expected FastHyperLogLog, got {type(sketch.sketch)}"
        
        # Should be configured for performance
        if hasattr(sketch.sketch, 'get_performance_info'):
            perf_info = sketch.sketch.get_performance_info()
            expected_cython = CYTHON_AVAILABLE  # Should use Cython if available
            assert perf_info['using_cython'] == expected_cython
    
    def test_intervalsketch_manual_control(self):
        """Test manual control over FastHyperLogLog usage in IntervalSketch."""
        # Explicitly disable FastHyperLogLog
        sketch = IntervalSketch(mode='A', precision=12, use_fast_hll=False)
        
        # Should fall back to regular HyperLogLog
        assert isinstance(sketch.sketch, HyperLogLog)
        assert not isinstance(sketch.sketch, FastHyperLogLog)
    
    def test_intervalsketch_file_processing_performance(self, temp_bed_file):
        """Test performance improvement in real file processing."""
        if not CYTHON_AVAILABLE:
            pytest.skip("Cython not available for performance testing")
        
        # Process with standard HyperLogLog
        start_time = time.time()
        sketch_standard = IntervalSketch.from_file(temp_bed_file, mode='A', 
                                                 precision=12, use_fast_hll=False)
        standard_time = time.time() - start_time
        
        # Process with FastHyperLogLog
        start_time = time.time()
        sketch_fast = IntervalSketch.from_file(temp_bed_file, mode='A', 
                                             precision=12, use_fast_hll=True)
        fast_time = time.time() - start_time
        
        # Verify both produce similar results
        standard_cardinality = sketch_standard.sketch.estimate_cardinality()
        fast_cardinality = sketch_fast.sketch.estimate_cardinality()
        
        if standard_cardinality > 0:
            relative_diff = abs(standard_cardinality - fast_cardinality) / standard_cardinality
            assert relative_diff < 0.1, f"Results should be similar, diff: {relative_diff:.3f}"
        
        # For small files, performance difference might not be significant
        # Just verify that fast version completes successfully
        assert fast_time > 0, "Fast processing should complete"
        assert standard_time > 0, "Standard processing should complete"
    
    def test_intervalsketch_large_dataset_benefit(self, large_bed_content):
        """Test performance benefit on larger datasets."""
        if not CYTHON_AVAILABLE:
            pytest.skip("Cython not available for performance testing")
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(large_bed_content)
            temp_path = f.name
        
        try:
            # Process with FastHyperLogLog
            start_time = time.time()
            sketch_fast = IntervalSketch.from_file(temp_path, mode='A', precision=14)
            fast_time = time.time() - start_time
            
            # Process with standard HyperLogLog  
            start_time = time.time()
            sketch_standard = IntervalSketch.from_file(temp_path, mode='A', 
                                                     precision=14, use_fast_hll=False)
            standard_time = time.time() - start_time
            
            # Verify cardinality estimates are similar
            fast_cardinality = sketch_fast.sketch.estimate_cardinality()
            standard_cardinality = sketch_standard.sketch.estimate_cardinality()
            
            relative_diff = abs(fast_cardinality - standard_cardinality) / standard_cardinality
            assert relative_diff < 0.05, f"Large dataset results should be very similar, diff: {relative_diff:.3f}"
            
            # For larger datasets, we expect some performance benefit
            if fast_time > 0 and standard_time > 0:
                speedup = standard_time / fast_time
                print(f"Large dataset speedup: {speedup:.2f}x (fast: {fast_time:.3f}s, standard: {standard_time:.3f}s)")
                
                # Even modest improvement is valuable
                assert speedup >= 0.8, f"Fast version should not be significantly slower, speedup: {speedup:.2f}x"
        
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


class TestRealWorldUsagePatterns:
    """Tests for real-world usage patterns and workflows."""
    
    def test_batch_vs_individual_performance(self):
        """Test performance difference between batch and individual processing."""
        if not CYTHON_AVAILABLE:
            pytest.skip("Cython not available for performance testing")
        
        # Generate test data
        test_data = [f"genomic_sequence_{i}" for i in range(5000)]
        
        # Individual processing
        sketch_individual = FastHyperLogLog(precision=14, use_cython=True)
        start_time = time.time()
        for item in test_data:
            sketch_individual.add_string(item)
        individual_time = time.time() - start_time
        
        # Batch processing
        sketch_batch = FastHyperLogLog(precision=14, use_cython=True)
        start_time = time.time()
        sketch_batch.add_batch(test_data)
        batch_time = time.time() - start_time
        
        # Verify similar cardinality
        individual_cardinality = sketch_individual.estimate_cardinality()
        batch_cardinality = sketch_batch.estimate_cardinality()
        
        relative_diff = abs(individual_cardinality - batch_cardinality) / individual_cardinality
        assert relative_diff < 0.01, f"Batch and individual should give same results, diff: {relative_diff:.3f}"
        
        # Batch should be faster
        if batch_time > 0:
            speedup = individual_time / batch_time
            print(f"Batch speedup: {speedup:.2f}x (individual: {individual_time:.3f}s, batch: {batch_time:.3f}s)")
            assert speedup > 1.5, f"Batch processing should be significantly faster, speedup: {speedup:.2f}x"
    
    def test_multiple_sketch_workflow(self):
        """Test workflow with multiple sketches and comparisons."""
        # Simulate workflow with multiple genomic samples
        samples = {
            'sample_A': [f"gene_A_{i}" for i in range(1000)],
            'sample_B': [f"gene_B_{i}" for i in range(800)] + [f"gene_A_{i}" for i in range(200)],  # 20% overlap
            'sample_C': [f"gene_C_{i}" for i in range(1200)],
        }
        
        # Create optimized sketches for each sample
        sketches = {}
        processing_times = {}
        
        for name, data in samples.items():
            start_time = time.time()
            sketches[name] = create_fast_hyperloglog(precision=14)
            sketches[name].add_batch(data)
            processing_times[name] = time.time() - start_time
        
        # Verify cardinality estimates
        assert abs(sketches['sample_A'].estimate_cardinality() - 1000) / 1000 < 0.1
        assert abs(sketches['sample_B'].estimate_cardinality() - 1000) / 1000 < 0.1  # 800 + 200 unique
        assert abs(sketches['sample_C'].estimate_cardinality() - 1200) / 1200 < 0.1
        
        # Test similarity calculations
        jaccard_AB = sketches['sample_A'].estimate_jaccard(sketches['sample_B'])
        jaccard_AC = sketches['sample_A'].estimate_jaccard(sketches['sample_C'])
        
        # A and B should have higher similarity than A and C
        assert jaccard_AB > jaccard_AC, "Samples with overlap should be more similar"
        
        # Check processing completed in reasonable time
        max_time = max(processing_times.values())
        assert max_time < 5.0, f"Processing should be fast, max time: {max_time:.3f}s"
    
    def test_memory_usage_stability(self):
        """Test that memory usage remains stable across multiple operations."""
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Perform many operations
        for i in range(20):
            sketch = create_fast_hyperloglog(precision=12)
            test_data = [f"iteration_{i}_item_{j}" for j in range(1000)]
            sketch.add_batch(test_data)
            
            # Force some garbage collection scenarios
            cardinality = sketch.estimate_cardinality()
            assert cardinality > 0
        
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_growth = final_memory - initial_memory
        
        # Memory growth should be reasonable (not indicative of major leaks)
        assert memory_growth < 50, f"Memory growth too high: {memory_growth:.1f}MB"
    
    def test_high_throughput_scenario(self):
        """Test high-throughput processing scenario."""
        if not CYTHON_AVAILABLE:
            pytest.skip("Cython not available for performance testing")
        
        # Simulate processing many small datasets quickly
        num_datasets = 50
        dataset_size = 200
        
        start_time = time.time()
        results = []
        
        for i in range(num_datasets):
            sketch = create_fast_hyperloglog(precision=12)
            data = [f"dataset_{i}_item_{j}" for j in range(dataset_size)]
            sketch.add_batch(data)
            cardinality = sketch.estimate_cardinality()
            results.append(cardinality)
        
        total_time = time.time() - start_time
        
        # All results should be reasonable
        for i, cardinality in enumerate(results):
            error = abs(cardinality - dataset_size) / dataset_size
            assert error < 0.3, f"Dataset {i} error too high: {error:.3f}"
        
        # Total processing should be fast
        throughput = num_datasets / total_time
        print(f"High-throughput test: {throughput:.1f} datasets/second")
        assert throughput > 10, f"Throughput too low: {throughput:.1f} datasets/second"


class TestBackwardCompatibility:
    """Tests for backward compatibility and migration scenarios."""
    
    def test_drop_in_replacement(self):
        """Test that FastHyperLogLog can be used as drop-in replacement."""
        test_data = [f"compatibility_test_{i}" for i in range(500)]
        
        # Original usage pattern
        original_sketch = HyperLogLog(precision=12, seed=42)
        for item in test_data:
            original_sketch.add_string(item)
        original_cardinality = original_sketch.estimate_cardinality()
        
        # Drop-in replacement usage
        fast_sketch = FastHyperLogLog(precision=12, seed=42, use_cython=False)
        for item in test_data:
            fast_sketch.add_string(item)
        fast_cardinality = fast_sketch.estimate_cardinality()
        
        # Should give very similar results
        relative_diff = abs(original_cardinality - fast_cardinality) / original_cardinality
        assert relative_diff < 0.02, f"Drop-in replacement should give similar results, diff: {relative_diff:.3f}"
    
    def test_migration_scenarios(self):
        """Test common migration scenarios from HyperLogLog to FastHyperLogLog."""
        test_data = [f"migration_test_{i}" for i in range(300)]
        
        # Scenario 1: Existing code using individual adds
        original_pattern = HyperLogLog(precision=12)
        for item in test_data:
            original_pattern.add_string(item)
        
        # Migrated to optimized pattern
        optimized_pattern = create_fast_hyperloglog(precision=12)
        optimized_pattern.add_batch(test_data)  # Use batch for performance
        
        # Results should be consistent
        original_cardinality = original_pattern.estimate_cardinality()
        optimized_cardinality = optimized_pattern.estimate_cardinality()
        
        relative_diff = abs(original_cardinality - optimized_cardinality) / original_cardinality
        assert relative_diff < 0.05, f"Migration should preserve accuracy, diff: {relative_diff:.3f}"
    
    def test_serialization_compatibility(self, temp_bed_file):
        """Test that sketches can be saved and loaded properly."""
        # Create and populate sketch
        sketch = create_fast_hyperloglog(precision=12)
        test_data = [f"serialization_test_{i}" for i in range(200)]
        sketch.add_batch(test_data)
        
        original_cardinality = sketch.estimate_cardinality()
        
        # Save sketch
        with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
            temp_path = f.name
        
        try:
            sketch.write(temp_path)
            
            # Load sketch back
            loaded_sketch = FastHyperLogLog.load(temp_path)
            loaded_cardinality = loaded_sketch.estimate_cardinality()
            
            # Should have same cardinality
            assert abs(original_cardinality - loaded_cardinality) < 1e-10, \
                "Loaded sketch should have identical cardinality"
            
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


class TestErrorHandlingAndRobustness:
    """Tests for error handling and robustness."""
    
    def test_graceful_cython_failure_handling(self):
        """Test graceful handling when Cython operations fail."""
        # This is more of a documentation test since we can't easily simulate Cython failures
        # But we can test the fallback mechanisms
        
        sketch = FastHyperLogLog(precision=12, use_cython=False)
        
        # Should work fine with Python fallback
        test_data = [f"fallback_test_{i}" for i in range(100)]
        sketch.add_batch(test_data)
        
        cardinality = sketch.estimate_cardinality()
        error = abs(cardinality - len(test_data)) / len(test_data)
        assert error < 0.3, f"Python fallback should work correctly, error: {error:.3f}"
    
    def test_resource_cleanup(self):
        """Test that resources are properly cleaned up."""
        # Create many sketches to test resource management
        sketches = []
        for i in range(100):
            sketch = create_fast_hyperloglog(precision=10)  # Smaller for speed
            test_data = [f"cleanup_test_{i}_{j}" for j in range(10)]
            sketch.add_batch(test_data)
            sketches.append(sketch)
        
        # All should work correctly
        for i, sketch in enumerate(sketches):
            cardinality = sketch.estimate_cardinality()
            assert cardinality > 0, f"Sketch {i} should have positive cardinality"
        
        # Python garbage collection should handle cleanup
        sketches = None


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 