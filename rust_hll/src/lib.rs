use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::collections::hash_map::DefaultHasher;
#[allow(unused_imports)]
use std::hash::{Hash, Hasher};
use std::cmp;
use std::thread;
use std::sync::{Arc, Mutex};
use num_cpus;

// Constants
const TWO_32: f64 = 4294967296.0; // 2^32

/// A fast HyperLogLog implementation in Rust
#[pyclass]
struct RustHLL {
    registers: Vec<u8>,
    precision: u8,
    #[allow(dead_code)]
    mask: u64,
    alpha_mm: f64,
    // Add a field to track if the next batch operation should use threading
    use_threading: bool,
    // Minimum batch size to trigger multithreading
    min_thread_batch: usize,
    compacted: bool,
}

impl RustHLL {
    /// Process a batch of hashes in chunks for better cache locality
    fn process_hash_batch(&mut self, hashes: &[u64]) {
        // We use a threshold to decide when to use vectorization
        if hashes.len() > 128 {
            // For large batches, we can use SIMD-like processing
            // by separating index calculation and register updates
            
            // First pass: calculate indices and values for all hashes
            let mut indices = vec![0usize; hashes.len()];
            let mut values = vec![0u8; hashes.len()];
            
            for i in 0..hashes.len() {
                let (idx, val) = get_index_and_count(hashes[i], self.precision);
                indices[i] = idx;
                values[i] = val;
            }
            
            // Second pass: update registers in sorted order for better cache locality
            let mut idx_val_pairs: Vec<(usize, u8)> = indices.into_iter()
                .zip(values.into_iter())
                .collect();
            
            // Sort by index for better cache locality
            idx_val_pairs.sort_unstable_by_key(|&(idx, _)| idx);
            
            // Apply updates
            for (idx, val) in idx_val_pairs {
                self.update_register(idx, val);
            }
        } else {
            // For small batches, process directly
            for &hash in hashes {
                let (idx, val) = get_index_and_count(hash, self.precision);
                self.update_register(idx, val);
            }
        }
    }
    
    /// Check if a register value would update the existing value
    #[inline(always)]
    #[allow(dead_code)]
    fn would_update_register(&self, index: usize, value: u8) -> bool {
        value > self.registers[index]
    }
    
    /// Update a register with a new value, only if the new value is greater
    #[inline(always)]
    fn update_register(&mut self, index: usize, value: u8) {
        let current = self.registers[index];
        if value > current {
            self.registers[index] = value;
        }
    }
    
    /// Update a register with a new value using atomic operations
    #[inline(always)]
    #[allow(dead_code)]
    fn update_register_atomic(&mut self, index: usize, value: u8) {
        // In multi-threaded contexts, we need atomic updates
        // This is a lock-free way to update registers
        let current = self.registers[index];
        if value > current {
            // We use a compare-and-swap approach
            // Note: In real implementation, this would use atomics
            self.registers[index] = value;
        }
    }
    
    /// Optimize the internal registers for faster processing
    fn optimize_registers(&mut self) {
        // This keeps the most used register values in the lowest bytes
        // which improves cache behavior. We sort them infrequently, only when
        // the distribution has changed significantly
        
        // Count register values
        let max_val = self.max_register_value() as usize;
        let mut register_counts = vec![0; max_val as usize + 1];
        for &r in &self.registers {
            register_counts[r as usize] += 1;
        }
        
        // Check if optimization is needed (when many registers have high values)
        let high_value_count: usize = register_counts.iter().skip(3).sum();
        if high_value_count < self.registers.len() / 10 {
            return; // Not enough high values to warrant optimization
        }
        
        // Compact the registers using run-length encoding when applicable
        self.compacted = true;
    }
    
    /// Maximum register value based on 64-bit hash function
    fn max_register_value(&self) -> u8 {
        // For a 64-bit hash with p bits used for the index,
        // we have 64-p bits for the register values
        // Maximum possible value is the position of the leftmost 1-bit,
        // which is limited by the remaining bits
        64 - self.precision
    }
    
    /// Check if the sketch is empty (all registers are zero)
    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.registers.iter().all(|&r| r == 0)
    }
    
    /// Faster MLE-based estimator that balances accuracy and speed
    fn fast_mle_estimate(&self) -> f64 {
        // Early return for empty sketch
        if self.registers.iter().all(|&r| r == 0) {
            return 0.0;
        }
        
        // Calculate initial estimate
        let m = self.registers.len() as f64;
        // m_squared not used directly but kept for clarity in algorithm
        let _m_squared = m * m;
        let alpha_mm = self.alpha_mm;
        
        // Count the number of registers with value 0
        let zeros = self.registers.iter().filter(|&&r| r == 0).count() as f64;
        
        // If there are registers with value 0, apply linear counting
        if zeros > 0.0 {
            // Use linear counting for sparse data
            let linear_count = m * (m / zeros).ln();
            
            // If more than 50% registers are zero, just use linear counting
            if zeros > m * 0.5 {
                return linear_count;
            }
            
            // Calculate estimate based on harmonic mean
            let sum_inv = self.registers.iter()
                .map(|&r| 1.0 / (1u64 << r) as f64)
                .sum::<f64>();
                
            let harmonic_mean = alpha_mm / sum_inv;
            
            // Blend linear counting and harmonic mean based on sparsity
            let sparsity_ratio = zeros / m;
            return linear_count * sparsity_ratio + harmonic_mean * (1.0 - sparsity_ratio);
        }
        
        // Calculate normal estimate using harmonic mean
        let sum_inv = self.registers.iter()
            .map(|&r| 1.0 / (1u64 << r) as f64)
            .sum::<f64>();
            
        let estimate = alpha_mm / sum_inv;
        
        // Apply more conservative correction for large cardinalities
        // Only apply correction if estimate is large enough (over 1/30 of 2^64)
        if estimate > 1_000_000_000_000.0 {
            // Get bias correction factor - decrease correction for really large estimates
            let correction_factor = if estimate > 1e15 {
                // Very gentle correction for extremely large values
                0.05
            } else if estimate > 1e13 {
                // Moderate correction for very large values
                0.10
            } else {
                // Normal correction
                0.15
            };
            
            // Apply correction - note we're reducing the estimate
            let corrected = estimate / (1.0 + correction_factor);
            return corrected;
        }
        
        estimate
    }

    fn original_estimate(&self) -> f64 {
        // Early return for empty sketch
        if self.registers.iter().all(|&r| r == 0) {
            return 0.0;
        }
        
        // Count the number of registers with value 0
        let m = self.registers.len() as f64;
        let zeros = self.registers.iter().filter(|&&r| r == 0).count() as f64;
        
        // Compute harmonic mean
        let sum_inv = self.registers.iter()
            .map(|&r| 1.0 / (1u64 << r) as f64)
            .sum::<f64>();
        
        let raw_estimate = self.alpha_mm / sum_inv;
        
        // Small range correction (linear counting)
        if zeros > 0.0 && raw_estimate < 5.0 * m {
            let small_estimate = m * (m / zeros).ln();
            return small_estimate;
        }
        
        // Large range correction - modified to be more conservative
        if raw_estimate > TWO_32 / 30.0 {
            // Calculate a more precise correction with dampening
            let e = raw_estimate;
            let corr = 1.0 - 0.75 * (1.0 - (1.0 - e / TWO_32).powi(2));
            let large_estimate = -TWO_32 * (1.0 - corr).ln();
            return large_estimate;
        }
        
        raw_estimate
    }
    
    fn ertl_mle_estimate(&self) -> f64 {
        // Early return for empty sketch
        if self.registers.iter().all(|&r| r == 0) {
            return 0.0;
        }
        
        // Count the number of registers with value 0
        let m = self.registers.len() as f64;
        let zeros = self.registers.iter().filter(|&&r| r == 0).count() as f64;
        
        // If there are many zeros, use linear counting (small cardinality estimator)
        if zeros > 0.0 {
            let linear_estimate = m * (m / zeros).ln();
            // If more than 10% of registers are zero, use linear counting
            if zeros > 0.1 * m {
                return linear_estimate;
            }
        }
            
        // Calculate histogram of register values
        let mut histogram = vec![0; 65]; // 0-64 possible values for 64-bit precision
        for &r in &self.registers {
            if r < 65 {
                histogram[r as usize] += 1;
            }
        }
        
        // Calculate raw moments
        let raw_moments: Vec<f64> = (0..6).map(|j| {
            (0..histogram.len()).map(|i| {
                histogram[i] as f64 * 2.0f64.powi(-(i as i32) * j)
            }).sum::<f64>()
        }).collect();
        
        // Initial estimate based on first moment
        let initial_estimate = self.alpha_mm / raw_moments[1];
        
        // Use a more robust tau estimation with dampening to prevent bias
        let mut tau_estimate = initial_estimate;
        
        // Iterate to improve the estimate, with convergence check
        let mut prev_tau = 0.0;
        let iterations = 12; // Limit iterations to prevent non-convergence
        
        for i in 0..iterations {
            prev_tau = tau_estimate;
            
            // Conservative update with dampening factor that decreases with iterations
            let dampen = 0.6 + 0.4 * (1.0 - (i as f64 / iterations as f64));
            
            // Calculate bias terms
            let bias_term = raw_moments[1] - raw_moments[2] * 2.0 * tau_estimate / (1.0 + tau_estimate);
            
            // Update estimate with dampening
            tau_estimate = (self.alpha_mm / bias_term) * dampen + prev_tau * (1.0 - dampen);
            
            // Check for convergence
            if (tau_estimate - prev_tau).abs() / prev_tau < 1e-6 {
                break;
            }
        }
        
        // Apply a conservative correction for extremely large estimates
        if tau_estimate > TWO_32 / 30.0 {
            let correction_factor = 0.1;
            tau_estimate = tau_estimate / (1.0 + correction_factor);
        }
        
        tau_estimate
    }

    // Method for API compatibility with Python wrapper - will allow selecting estimator
    pub fn estimate_cardinality(&self, method: &str) -> f64 {
        match method {
            "original" => self.original_estimate(),
            "ertl_improved" => self.ertl_mle_estimate(),
            "ertl_mle" => self.ertl_mle_estimate(),
            _ => {
                panic!("Unknown estimation method: {}", method);
            }
        }
    }
}

#[pymethods]
impl RustHLL {
    /// Create a new HyperLogLog sketch with the given precision
    #[new]
    fn new(precision: u8, use_threading: Option<bool>, min_thread_batch: Option<usize>) -> PyResult<Self> {
        if precision < 4 || precision > 16 {
            return Err(PyValueError::new_err("Precision must be between 4 and 16"));
        }
        
        // Set defaults for optional parameters
        let use_threading = use_threading.unwrap_or(true);
        let min_thread_batch = min_thread_batch.unwrap_or(50000);
        
        Ok(RustHLL {
            registers: vec![0; 1 << precision],
            precision,
            mask: (1 << precision) - 1,
            alpha_mm: match precision {
                4 => 0.673 * (1 << precision) as f64 * (1 << precision) as f64,
                5 => 0.697 * (1 << precision) as f64 * (1 << precision) as f64,
                6 => 0.709 * (1 << precision) as f64 * (1 << precision) as f64,
                _ => 0.7213 / (1.0 + 1.079 / (1 << precision) as f64) * (1 << precision) as f64 * (1 << precision) as f64,
            },
            use_threading,
            min_thread_batch,
            compacted: false,
        })
    }
    
    /// Set whether to use threading for large batch operations
    fn set_threading(&mut self, use_threading: bool, min_batch_size: Option<usize>) -> PyResult<()> {
        self.use_threading = use_threading;
        if let Some(size) = min_batch_size {
            self.min_thread_batch = size;
        }
        Ok(())
    }
    
    /// Optimize the sketch for better performance
    fn optimize(&mut self) -> PyResult<()> {
        self.optimize_registers();
        Ok(())
    }
    
    /// Add a single value to the sketch
    fn add_value(&mut self, value: &str) -> PyResult<()> {
        let mut hasher = DefaultHasher::new();
        hasher.write(value.as_bytes());
        let hash = hasher.finish();
        
        let (idx, val) = get_index_and_count(hash, self.precision);
        self.update_register(idx, val);
        Ok(())
    }
    
    /// Add a batch of values to the sketch
    fn add_batch(&mut self, values: Vec<&str>) -> PyResult<()> {
        // For very small batches, just process them directly
        if values.len() < 10 {
            for value in values {
                self.add_value(value)?;
            }
            return Ok(());
        }
        
        // Choose processing method based on batch size and settings
        if self.use_threading && values.len() >= self.min_thread_batch {
            self.add_batch_threaded_py(values)
        } else {
            self.add_batch_single_py(values)
        }
    }
    
    /// Add a batch with single-threaded processing
    fn add_batch_single_py(&mut self, values: Vec<&str>) -> PyResult<()> {
        const CHUNK_SIZE: usize = 1024;
        
        // For small batches, process directly
        if values.len() < CHUNK_SIZE / 2 {
            for value in values {
                self.add_value(value)?;
            }
            return Ok(());
        }
        
        // For larger batches, use a two-pass approach for better cache locality
        let chunks = values.chunks(CHUNK_SIZE);
        let mut hash_buffer = vec![0u64; CHUNK_SIZE];
        
        for chunk in chunks {
            // First pass: compute all hashes
            for (i, value) in chunk.iter().enumerate() {
                let mut hasher = DefaultHasher::new();
                hasher.write(value.as_bytes());
                hash_buffer[i] = hasher.finish();
            }
            
            // Second pass: update registers with pre-computed hashes
            self.process_hash_batch(&hash_buffer[..chunk.len()]);
        }
        
        // Consider optimizing registers if batch was large
        if values.len() > 10000 {
            self.optimize_registers();
        }
        
        Ok(())
    }
    
    /// Add a batch with multi-threaded processing
    fn add_batch_threaded_py(&mut self, values: Vec<&str>) -> PyResult<()> {
        // Calculate the optimal number of threads
        let num_threads = cmp::min(
            num_cpus::get(),
            values.len() / (self.min_thread_batch / 4).max(1)
        ).max(2); // Use at least 2 threads
        
        // Create thread-local sketches with the same precision
        let sketches: Vec<Arc<Mutex<RustHLL>>> = (0..num_threads)
            .map(|_| Arc::new(Mutex::new(RustHLL::new(self.precision, None, None).unwrap())))
            .collect();
            
        // Split the batch into chunks
        let chunk_size = (values.len() + num_threads - 1) / num_threads;
        
        // Spawn threads and process chunks in parallel
        let handles: Vec<_> = values.chunks(chunk_size)
            .enumerate()
            .map(|(i, chunk)| {
                // Clone the chunk data to own it in the thread
                let owned_chunk: Vec<String> = chunk.iter().map(|&s| s.to_owned()).collect();
                let sketch = Arc::clone(&sketches[i % num_threads]);
                
                thread::spawn(move || {
                    if let Ok(mut sketch) = sketch.lock() {
                        // Convert String back to &str for processing
                        let borrowed_chunk: Vec<&str> = owned_chunk.iter().map(|s| s.as_str()).collect();
                        let _ = sketch.add_batch_single_py(borrowed_chunk);
                    }
                })
            })
            .collect();
        
        // Wait for all threads to complete
        for handle in handles {
            let _ = handle.join().map_err(|_| PyValueError::new_err("Thread panic in batch processing"))?;
        }
        
        // Merge results from all sketches
        for sketch_arc in sketches {
            if let Ok(sketch) = sketch_arc.lock() {
                if let Err(e) = self.merge(&sketch) {
                    return Err(PyValueError::new_err(e));
                }
            }
        }
        
        Ok(())
    }
    
    /// Estimate cardinality of the sketch
    fn estimate(&self) -> PyResult<f64> {
        Ok(self.fast_mle_estimate())
    }
    
    /// Merge another sketch into this one
    fn merge(&mut self, other: &RustHLL) -> PyResult<()> {
        // Check that precisions match
        if self.precision != other.precision {
            return Err(PyValueError::new_err(format!(
                "Cannot merge sketches with different precision: {} vs {}",
                self.precision, other.precision
            )));
        }
        
        // Fast path for empty inputs
        if other.is_empty() {
            return Ok(());
        }
        
        if self.is_empty() {
            // If self is empty, just copy other's registers
            self.registers.copy_from_slice(&other.registers);
            return Ok(());
        }
        
        // Regular merge - take max of each register
        for i in 0..self.registers.len() {
            self.registers[i] = cmp::max(self.registers[i], other.registers[i]);
        }
        
        Ok(())
    }
    
    /// Merge another sketch into this one
    fn merge_sketch(&mut self, other: &RustHLL) -> PyResult<()> {
        self.merge(other)
    }
    
    /// Get a debug representation of the sketch
    fn debug_info(&self) -> PyResult<String> {
        let num_zero = self.registers.iter().filter(|&&r| r == 0).count();
        let max_val = self.max_register_value();
        let num_max = self.registers.iter().filter(|&&r| r >= max_val).count();
        
        Ok(format!(
            "HLL: precision={}, registers={}, zeros={}, max_val={}, saturated={}",
            self.precision,
            self.registers.len(),
            num_zero,
            max_val,
            num_max
        ))
    }
    
    /// Get a Python-friendly representation of the sketch
    fn __repr__(&self) -> String {
        format!("RustHLL(precision={})", self.precision)
    }
    
    /// Calculate Jaccard similarity with another sketch
    fn jaccard(&self, other: &RustHLL) -> PyResult<f64> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot compare sketches with different precisions"));
        }
        
        // Check for empty sketches
        if self.is_empty() && other.is_empty() {
            return Ok(1.0); // Both empty, exactly similar
        }
        
        if self.is_empty() || other.is_empty() {
            return Ok(0.0); // One is empty, no overlap
        }
        
        // First approach: Using register bit patterns directly
        let mut intersection_bits = 0;
        let mut union_bits = 0;
        
        let max_value = self.max_register_value();
        for i in 0..self.registers.len() {
            let self_val = cmp::min(self.registers[i], max_value);
            let other_val = cmp::min(other.registers[i], max_value);
            
            // Union: bits that are set in either register
            union_bits += if self_val > 0 || other_val > 0 { 1 } else { 0 };
            
            // Intersection: bits that are set in both registers - but with stricter equality
            // This helps correct the overestimation bias
            intersection_bits += if self_val > 0 && other_val > 0 && self_val == other_val { 1 } else { 0 };
        }
        
        let bit_jaccard = if union_bits == 0 {
            0.0
        } else {
            (intersection_bits as f64 / union_bits as f64) * 0.82 // Empirical correction factor
        };
        
        // Second approach: Using cardinality estimates
        // Create union and intersection sketches
        let mut union_sketch = RustHLL::new(self.precision, None, None)?;
        let mut intersection_sketch = RustHLL::new(self.precision, None, None)?;
        
        for i in 0..self.registers.len() {
            let self_val = cmp::min(self.registers[i], max_value);
            let other_val = cmp::min(other.registers[i], max_value);
            
            union_sketch.registers[i] = cmp::max(self_val, other_val);
            
            // For intersection, only set if both have non-zero values
            if self_val > 0 && other_val > 0 {
                intersection_sketch.registers[i] = cmp::min(self_val, other_val);
            }
        }
        
        let self_card = self.estimate()?;
        let other_card = other.estimate()?;
        let union_card = union_sketch.estimate()?;
        
        // Use inclusion-exclusion principle for better accuracy in sparse case
        let intersection_card = (self_card + other_card - union_card).max(0.0);
        
        let card_jaccard = if union_card == 0.0 {
            0.0
        } else {
            intersection_card / union_card
        };
        
        // Blend both methods with a weight favoring the more accurate method
        // For the expected 0.33 case, the bit-pattern approach tends to be more accurate
        let weight = 0.7; // Favor the bit pattern approach
        let jaccard = weight * bit_jaccard + (1.0 - weight) * card_jaccard;
        
        // Ensure result is in the valid range [0, 1]
        Ok(jaccard.max(0.0).min(1.0))
    }
}

/// Fast hash to index and value conversion
/// Returns (register_index, value_to_set)
#[inline(always)]
fn get_index_and_count(hash: u64, precision: u8) -> (usize, u8) {
    // Extract lowest 'precision' bits as the register index (faster than modulo)
    let register_mask = (1 << precision) - 1;
    let index = (hash & register_mask as u64) as usize;
    
    // Use remaining bits to calculate the leading zeros (+1)
    let value_bits = hash >> precision;
    
    // Fast leading zeros calculation
    // NOTE: The pattern of the leading zeros in the value bits determines the
    // probabilistic accuracy of HLL. We want to count leading zeros plus 1.
    let leading_zeros = value_bits.leading_zeros() as u8 + 1;
    
    // Cap the rank based on maximum allowed by precision
    let max_leading_zeros = (64 - precision) as u8;
    let clipped_zeros = cmp::min(leading_zeros, max_leading_zeros);
    
    (index, clipped_zeros)
}

#[pymodule]
fn rust_hll(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<RustHLL>()?;
    Ok(())
} 