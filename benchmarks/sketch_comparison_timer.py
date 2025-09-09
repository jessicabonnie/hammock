#!/usr/bin/env python3
"""
Sketch Comparison Timer

This module provides utilities to measure the compute time for sketch comparison
operations separately from sketch creation. This is useful for benchmarking
the performance of different comparison algorithms without the overhead of
sketch construction.

Author: AI Assistant
Date: 2024
"""

import time
import psutil
from typing import Dict, List, Tuple, Any, Union, Optional
from contextlib import contextmanager
import statistics
from dataclasses import dataclass
from hammock.lib.abstractsketch import AbstractSketch


@dataclass
class TimingResult:
    """Container for timing results."""
    wall_time: float
    cpu_time: float
    memory_delta: float
    num_comparisons: int
    comparison_type: str
    sketch_type: str


class SketchComparisonTimer:
    """
    Utility class for measuring sketch comparison performance.
    
    This class provides methods to time sketch comparison operations
    separately from sketch creation, allowing for more accurate
    benchmarking of comparison algorithms.
    """
    
    def __init__(self, enable_memory_monitoring: bool = True):
        """
        Initialize the timer.
        
        Args:
            enable_memory_monitoring: Whether to track memory usage during comparisons
        """
        self.enable_memory_monitoring = enable_memory_monitoring
        self.process = psutil.Process() if enable_memory_monitoring else None
    
    @contextmanager
    def time_comparison(self, comparison_type: str = "jaccard"):
        """
        Context manager for timing a single comparison operation.
        
        Args:
            comparison_type: Type of comparison being performed (e.g., "jaccard", "intersection")
            
        Yields:
            TimingResult: Object containing timing information
        """
        # Get initial memory usage
        memory_start = self._get_memory_usage()
        
        # Start timing
        wall_start = time.perf_counter()
        cpu_start = time.process_time()
        
        # Yield control to the comparison operation
        result = TimingResult(
            wall_time=0.0,
            cpu_time=0.0,
            memory_delta=0.0,
            num_comparisons=1,
            comparison_type=comparison_type,
            sketch_type="unknown"
        )
        
        try:
            yield result
        finally:
            # End timing
            wall_end = time.perf_counter()
            cpu_end = time.process_time()
            memory_end = self._get_memory_usage()
            
            # Calculate times
            result.wall_time = wall_end - wall_start
            result.cpu_time = cpu_end - cpu_start
            result.memory_delta = memory_end - memory_start
    
    def time_jaccard_comparison(self, sketch1: AbstractSketch, sketch2: AbstractSketch, 
                              num_runs: int = 1) -> List[TimingResult]:
        """
        Time Jaccard similarity comparison between two sketches.
        
        Args:
            sketch1: First sketch for comparison
            sketch2: Second sketch for comparison
            num_runs: Number of times to run the comparison for averaging
            
        Returns:
            List of TimingResult objects, one for each run
        """
        results = []
        
        for run in range(num_runs):
            with self.time_comparison("jaccard") as timer:
                # Perform the comparison
                jaccard_value = sketch1.estimate_jaccard(sketch2)
                timer.sketch_type = f"{type(sketch1).__name__}_vs_{type(sketch2).__name__}"
                timer.jaccard_value = jaccard_value  # Store the result for reference
            
            results.append(timer)
        
        return results
    
    def time_intersection_comparison(self, sketch1: AbstractSketch, sketch2: AbstractSketch,
                                   num_runs: int = 1) -> List[TimingResult]:
        """
        Time intersection estimation between two sketches.
        
        Args:
            sketch1: First sketch for comparison
            sketch2: Second sketch for comparison
            num_runs: Number of times to run the comparison for averaging
            
        Returns:
            List of TimingResult objects, one for each run
        """
        results = []
        
        for run in range(num_runs):
            with self.time_comparison("intersection") as timer:
                # Perform the intersection estimation
                intersection_value = sketch1.estimate_intersection(sketch2)
                timer.sketch_type = f"{type(sketch1).__name__}_vs_{type(sketch2).__name__}"
                timer.intersection_value = intersection_value  # Store the result for reference
            
            results.append(timer)
        
        return results
    
    def time_union_comparison(self, sketch1: AbstractSketch, sketch2: AbstractSketch,
                            num_runs: int = 1) -> List[TimingResult]:
        """
        Time union estimation between two sketches.
        
        Args:
            sketch1: First sketch for comparison
            sketch2: Second sketch for comparison
            num_runs: Number of times to run the comparison for averaging
            
        Returns:
            List of TimingResult objects, one for each run
        """
        results = []
        
        for run in range(num_runs):
            with self.time_comparison("union") as timer:
                # Perform the union estimation
                union_value = sketch1.estimate_union(sketch2)
                timer.sketch_type = f"{type(sketch1).__name__}_vs_{type(sketch2).__name__}"
                timer.union_value = union_value  # Store the result for reference
            
            results.append(timer)
        
        return results
    
    def time_similarity_values(self, sketch1: AbstractSketch, sketch2: AbstractSketch,
                             num_runs: int = 1) -> List[TimingResult]:
        """
        Time similarity_values method between two sketches.
        
        Args:
            sketch1: First sketch for comparison
            sketch2: Second sketch for comparison
            num_runs: Number of times to run the comparison for averaging
            
        Returns:
            List of TimingResult objects, one for each run
        """
        results = []
        
        for run in range(num_runs):
            with self.time_comparison("similarity_values") as timer:
                # Perform the similarity calculation
                similarity_dict = sketch1.similarity_values(sketch2)
                timer.sketch_type = f"{type(sketch1).__name__}_vs_{type(sketch2).__name__}"
                timer.similarity_dict = similarity_dict  # Store the result for reference
            
            results.append(timer)
        
        return results
    
    def benchmark_comparison_methods(self, sketch1: AbstractSketch, sketch2: AbstractSketch,
                                   num_runs: int = 5) -> Dict[str, Dict[str, float]]:
        """
        Benchmark all available comparison methods between two sketches.
        
        Args:
            sketch1: First sketch for comparison
            sketch2: Second sketch for comparison
            num_runs: Number of runs for each comparison method
            
        Returns:
            Dictionary containing timing statistics for each comparison method
        """
        results = {}
        
        # Test Jaccard comparison
        if hasattr(sketch1, 'estimate_jaccard'):
            jaccard_results = self.time_jaccard_comparison(sketch1, sketch2, num_runs)
            results['jaccard'] = self._calculate_statistics(jaccard_results)
        
        # Test intersection comparison
        if hasattr(sketch1, 'estimate_intersection'):
            intersection_results = self.time_intersection_comparison(sketch1, sketch2, num_runs)
            results['intersection'] = self._calculate_statistics(intersection_results)
        
        # Test union comparison
        if hasattr(sketch1, 'estimate_union'):
            union_results = self.time_union_comparison(sketch1, sketch2, num_runs)
            results['union'] = self._calculate_statistics(union_results)
        
        # Test similarity_values
        if hasattr(sketch1, 'similarity_values'):
            similarity_results = self.time_similarity_values(sketch1, sketch2, num_runs)
            results['similarity_values'] = self._calculate_statistics(similarity_results)
        
        return results
    
    def _calculate_statistics(self, results: List[TimingResult]) -> Dict[str, float]:
        """Calculate timing statistics from a list of TimingResult objects."""
        if not results:
            return {}
        
        wall_times = [r.wall_time for r in results]
        cpu_times = [r.cpu_time for r in results]
        memory_deltas = [r.memory_delta for r in results]
        
        return {
            'wall_time_mean': statistics.mean(wall_times),
            'wall_time_std': statistics.stdev(wall_times) if len(wall_times) > 1 else 0.0,
            'wall_time_min': min(wall_times),
            'wall_time_max': max(wall_times),
            'cpu_time_mean': statistics.mean(cpu_times),
            'cpu_time_std': statistics.stdev(cpu_times) if len(cpu_times) > 1 else 0.0,
            'cpu_time_min': min(cpu_times),
            'cpu_time_max': max(cpu_times),
            'memory_delta_mean': statistics.mean(memory_deltas),
            'memory_delta_std': statistics.stdev(memory_deltas) if len(memory_deltas) > 1 else 0.0,
            'num_runs': len(results)
        }
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        if not self.enable_memory_monitoring or not self.process:
            return 0.0
        
        try:
            memory_info = self.process.memory_info()
            return memory_info.rss / (1024 * 1024)  # Convert to MB
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return 0.0
    
    def print_timing_summary(self, results: Dict[str, Dict[str, float]], 
                           sketch1_type: str, sketch2_type: str):
        """
        Print a formatted summary of timing results.
        
        Args:
            results: Results from benchmark_comparison_methods
            sketch1_type: Type of first sketch
            sketch2_type: Type of second sketch
        """
        print(f"\n{'='*80}")
        print(f"SKETCH COMPARISON TIMING SUMMARY")
        print(f"Sketch Types: {sketch1_type} vs {sketch2_type}")
        print(f"{'='*80}")
        
        for method, stats in results.items():
            print(f"\n{method.upper()} COMPARISON:")
            print(f"  Wall Time: {stats['wall_time_mean']:.6f}s ± {stats['wall_time_std']:.6f}s")
            print(f"    Range: {stats['wall_time_min']:.6f}s - {stats['wall_time_max']:.6f}s")
            print(f"  CPU Time:  {stats['cpu_time_mean']:.6f}s ± {stats['cpu_time_std']:.6f}s")
            print(f"    Range: {stats['cpu_time_min']:.6f}s - {stats['cpu_time_max']:.6f}s")
            print(f"  Memory:    {stats['memory_delta_mean']:.2f}MB ± {stats['memory_delta_std']:.2f}MB")
            print(f"  Runs:      {stats['num_runs']}")


def compare_sketch_types(sketch_pairs: List[Tuple[AbstractSketch, AbstractSketch, str, str]],
                        num_runs: int = 5) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Compare timing performance across different sketch type pairs.
    
    Args:
        sketch_pairs: List of (sketch1, sketch2, name1, name2) tuples
        num_runs: Number of runs for each comparison
        
    Returns:
        Dictionary mapping sketch pair names to timing results
    """
    timer = SketchComparisonTimer()
    all_results = {}
    
    for sketch1, sketch2, name1, name2 in sketch_pairs:
        pair_name = f"{name1}_vs_{name2}"
        print(f"\nBenchmarking {pair_name}...")
        
        results = timer.benchmark_comparison_methods(sketch1, sketch2, num_runs)
        all_results[pair_name] = results
        
        timer.print_timing_summary(results, name1, name2)
    
    return all_results


if __name__ == "__main__":
    # Example usage
    print("Sketch Comparison Timer - Example Usage")
    print("This module provides utilities for timing sketch comparison operations.")
    print("Import this module and use SketchComparisonTimer class for benchmarking.")
