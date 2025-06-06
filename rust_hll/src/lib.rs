use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use ahash::{AHasher, RandomState};
use std::hash::{Hash, Hasher, BuildHasher};
use std::cmp;
use std::thread;
use std::sync::{Arc, Mutex};
use num_cpus;
use bincode;
use serde::{Serialize, Deserialize};
use std::cell::RefCell;
use std::sync::atomic::{AtomicU8, Ordering};

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
    // Track hash size (32 or 64 bits)
    hash_size: u8,
    hasher: RandomState,
    memory_limit: Option<usize>,  // Add memory limit field
    debug: bool,  // Add debug field
    seed: Option<u64>,  // Add seed field
}

// Internal implementation
impl RustHLL {
    /// Process a batch of hashes in chunks for better cache locality and memory efficiency
    fn process_hash_batch(&mut self, hashes: &[u64]) {
        if self.debug {
            println!("Processing batch of {} hashes", hashes.len());
        }
        
        // Process in smaller chunks to limit memory usage
        const CHUNK_SIZE: usize = 256;  // Increased back to original value
        
        // Pre-allocate hash vector with fixed capacity
        let mut hash_buffer = Vec::with_capacity(CHUNK_SIZE);
            
        for (i, chunk) in hashes.chunks(CHUNK_SIZE).enumerate() {
            if self.debug && i % 100 == 0 {
                println!("Processing chunk {}", i);
            }
            
            // Clear and reuse the buffer
            hash_buffer.clear();
            hash_buffer.extend_from_slice(chunk);
            
            // Process each chunk directly without creating temporary arrays
            for &hash in &hash_buffer {
                let (idx, val) = get_index_and_count(hash, self.precision, self.hash_size);
                if self.registers[idx] < val {
                    self.registers[idx] = val;
                }
            }
        }
        
        if self.debug {
            println!("Finished processing batch");
        }
    }

    fn would_update_register(&self, index: usize, value: u8) -> bool {
        self.registers[index] < value
    }

    fn update_register(&mut self, index: usize, value: u8) {
        if self.would_update_register(index, value) {
            self.registers[index] = value;
        }
    }

    #[allow(dead_code)]
    fn update_register_atomic(&mut self, index: usize, value: u8) {
        if self.would_update_register(index, value) {
            self.registers[index] = value;
        }
    }

    fn optimize_registers(&mut self) {
        if !self.compacted {
            let max_val = self.max_register_value();
            // More aggressive compaction strategy
            if max_val > 20 {  // Reduced from 30 to 20
                // Compact registers if max value is too large
                for r in &mut self.registers {
                    if *r > 20 {
                        *r = 20;
                    }
                }
                self.compacted = true;
            }
        }
    }

    fn max_register_value(&self) -> u8 {
        self.registers.iter().max().copied().unwrap_or(0)
    }

    fn fast_mle_estimate(&self) -> f64 {
        let sum: f64 = self.registers.iter()
            .map(|&r| 1.0 / (1 << r) as f64)
            .sum();
        
        let m = self.registers.len() as f64;
        let alpha_mm = self.alpha_mm;
        
        let estimate = alpha_mm * m * m / sum;
        
        // Apply bias correction for small cardinalities
        if estimate <= 2.5 * m {
            let zeros = self.registers.iter().filter(|&&r| r == 0).count() as f64;
            if zeros > 0.0 {
                m * (m / zeros).ln()
            } else {
                estimate
            }
        } else {
            estimate
        }
    }

    #[allow(dead_code)]
    fn original_estimate(&self) -> f64 {
        let sum: f64 = self.registers.iter()
            .map(|&r| 1.0 / (1 << r) as f64)
            .sum();
        
        let m = self.registers.len() as f64;
        let alpha_mm = self.alpha_mm;
        
        alpha_mm * m * m / sum
    }

    #[allow(dead_code)]
    fn ertl_mle_estimate(&self) -> f64 {
        let sum: f64 = self.registers.iter()
            .map(|&r| 1.0 / (1 << r) as f64)
            .sum();
        
        let m = self.registers.len() as f64;
        let alpha_mm = self.alpha_mm;
        
        let estimate = alpha_mm * m * m / sum;
        
        // Apply bias correction for small cardinalities
        if estimate <= 2.5 * m {
            let zeros = self.registers.iter().filter(|&&r| r == 0).count() as f64;
            if zeros > 0.0 {
                m * (m / zeros).ln()
            } else {
                estimate
            }
        } else {
            estimate
        }
    }
}

// Python interface implementation
#[pymethods]
impl RustHLL {
    #[new]
    #[pyo3(signature = (precision, use_threading=None, min_thread_batch=None, hash_size=None, memory_limit=None, debug=None, seed=None))]
    pub fn new(precision: u8, use_threading: Option<bool>, min_thread_batch: Option<usize>, hash_size: Option<u8>, memory_limit: Option<usize>, debug: Option<bool>, seed: Option<u64>) -> PyResult<Self> {
        let debug = debug.unwrap_or(false);
        
        if debug {
            println!("Creating new RustHLL with parameters:");
            println!("  - Precision: {}", precision);
            println!("  - Use threading: {:?}", use_threading);
            println!("  - Min thread batch: {:?}", min_thread_batch);
            println!("  - Hash size: {:?}", hash_size);
            println!("  - Memory limit: {:?}", memory_limit);
            println!("  - Debug: {}", debug);
            println!("  - Seed: {:?}", seed);
        }
        
        // Validate hash_size
        let hash_size = hash_size.unwrap_or(32);
        if hash_size != 32 && hash_size != 64 {
            return Err(PyValueError::new_err("hash_size must be 32 or 64"));
        }

        // Check if precision would exceed the hash length
        if precision >= hash_size {
            return Err(PyValueError::new_err(format!("Precision must be < {} ({}-bit hash)", hash_size, hash_size)));
        }

        let m = 1 << precision;
        
        // Check memory usage against limit
        let estimated_memory = m * std::mem::size_of::<u8>();
        
        // Ensure minimum memory limit of 1MB
        let min_memory = 1024 * 1024;  // 1MB
        
        // If memory limit is provided, ensure it's at least the minimum
        let memory_limit = memory_limit.map(|limit| {
            let adjusted_limit = std::cmp::max(limit, min_memory);
            if debug {
                println!("Rust memory limit adjustment:");
                println!("  - Original limit: {} bytes", limit);
                println!("  - Minimum memory: {} bytes", min_memory);
                println!("  - Adjusted limit: {} bytes", adjusted_limit);
            }
            adjusted_limit
        });
        
        if let Some(limit) = memory_limit {
            // Add 10% buffer to estimated memory for safety
            let estimated_memory_with_buffer = (estimated_memory as f64 * 1.1) as usize;
            
            if debug {
                println!("Rust memory check:");
                println!("  - Estimated memory: {} bytes", estimated_memory);
                println!("  - Estimated with buffer: {} bytes", estimated_memory_with_buffer);
            }
            
            if estimated_memory_with_buffer > limit {
                return Err(PyValueError::new_err(format!(
                    "Memory limit exceeded: {} bytes needed, {} bytes limit",
                    estimated_memory_with_buffer, limit
                )));
            }
        }
        
        // Calculate alpha_mm
        let alpha_mm = match hash_size {
            32 => match m {
                16384 => 0.673,  // 2^14
                32768 => 0.697,  // 2^15
                65536 => 0.709,  // 2^16
                _ => 0.7213 / (1.0 + 1.079 / m as f64),
            },
            64 => match m {
                16384 => 0.351,  // 2^14
                32768 => 0.532,  // 2^15
                65536 => 0.625,  // 2^16
                _ => 0.7213 / (1.0 + 1.079 / m as f64) * 0.85,  // Additional scaling for 64-bit
            },
            _ => return Err(PyValueError::new_err("hash_size must be 32 or 64")),
        };
        
        // Create hasher with seed if provided
        let hasher = if let Some(seed) = seed {
            RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3)
        } else {
            RandomState::new()
        };
        
        Ok(RustHLL {
            registers: vec![0; m],
            precision,
            mask: (1 << precision) - 1,
            alpha_mm,
            use_threading: use_threading.unwrap_or(true),
            min_thread_batch: min_thread_batch.unwrap_or(10000),
            compacted: false,
            hash_size,
            hasher,
            memory_limit,
            debug,
            seed,
        })
    }

    #[getter]
    pub fn hash_size(&self) -> u8 {
        self.hash_size
    }

    pub fn set_threading(&mut self, use_threading: bool, min_batch_size: Option<usize>) -> PyResult<()> {
        self.use_threading = use_threading;
        if let Some(size) = min_batch_size {
            self.min_thread_batch = size;
        }
        Ok(())
    }

    pub fn optimize(&mut self) -> PyResult<()> {
        self.optimize_registers();
        Ok(())
    }

    pub fn add_value(&mut self, value: &str) -> PyResult<()> {
        let mut hasher = self.hasher.build_hasher();
        value.hash(&mut hasher);
        let hash = if self.hash_size == 32 {
            hasher.finish() as u32 as u64
        } else {
            hasher.finish()
        };
        
        let (idx, val) = get_index_and_count(hash, self.precision, self.hash_size);
        self.update_register(idx, val);
        Ok(())
    }

    fn add_batch_single_py(&mut self, values: Vec<&str>) -> PyResult<()> {
        let mut hashes = Vec::with_capacity(values.len());
        
        // First pass: compute all hashes
        for value in values {
            let mut hasher = self.hasher.build_hasher();
            value.hash(&mut hasher);
            let hash = if self.hash_size == 32 {
                hasher.finish() as u32 as u64
            } else {
                hasher.finish()
            };
            hashes.push(hash);
        }
        
        // Second pass: process hashes in batches
        self.process_hash_batch(&hashes);
        Ok(())
    }

    fn add_batch_threaded_py(&mut self, values: Vec<String>) -> PyResult<()> {
        // Check memory limit before proceeding
        if let Some(limit) = self.memory_limit {
            let current_usage = self.registers.len() * std::mem::size_of::<u8>();
            let estimated_batch_usage = values.len() * std::mem::size_of::<String>();
            let total_estimated = current_usage + estimated_batch_usage;
            
            if self.debug {
                println!("Memory usage check:");
                println!("  - Current usage: {} bytes", current_usage);
                println!("  - Batch usage: {} bytes", estimated_batch_usage);
                println!("  - Total estimated: {} bytes", total_estimated);
                println!("  - Memory limit: {} bytes", limit);
            }
            
            if total_estimated > limit {
                return Err(PyValueError::new_err(format!(
                    "Memory limit would be exceeded. Current: {} bytes, Batch: {} bytes, Total: {} bytes, Limit: {} bytes",
                    current_usage, estimated_batch_usage, total_estimated, limit
                )));
            }
        }

        // For small batches, use single-threaded processing
        if values.len() < self.min_thread_batch {
            return self.add_batch_single_py(values.iter().map(|s| s.as_str()).collect());
        }

        // Process in chunks without threading
        const CHUNK_SIZE: usize = 256;  // Increased back to original value
        
        for chunk in values.chunks(CHUNK_SIZE) {
            let mut hashes = Vec::with_capacity(chunk.len());
            
            // Compute hashes for the chunk
            for value in chunk {
                let mut hasher = self.hasher.build_hasher();
                value.hash(&mut hasher);
                let hash = if self.hash_size == 32 {
                    hasher.finish() as u32 as u64
                } else {
                    hasher.finish()
                };
                hashes.push(hash);
            }
            
            // Process the chunk of hashes
            self.process_hash_batch(&hashes);
        }
        
        Ok(())
    }

    pub fn add_batch(&mut self, values: Vec<&str>) -> PyResult<()> {
        if values.len() >= self.min_thread_batch && self.use_threading {
            // Convert references to owned strings
            let owned_values: Vec<String> = values.iter().map(|s| s.to_string()).collect();
            self.add_batch_threaded_py(owned_values)
        } else {
            self.add_batch_single_py(values)
        }
    }

    pub fn estimate_cardinality(&self) -> PyResult<f64> {
        Ok(self.fast_mle_estimate())
    }

    pub fn merge(&mut self, other: &RustHLL) -> PyResult<()> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot merge sketches with different precisions"));
        }
        
        for i in 0..self.registers.len() {
            if self.registers[i] < other.registers[i] {
                self.registers[i] = other.registers[i];
            }
        }
        Ok(())
    }

    pub fn merge_sketch(&mut self, other: &RustHLL) -> PyResult<()> {
        self.merge(other)
    }

    pub fn debug_info(&self) -> PyResult<String> {
        Ok(format!(
            "RustHLL(precision={}, registers={}, max_value={})",
            self.precision,
            self.registers.len(),
            self.max_register_value()
        ))
    }

    pub fn __repr__(&self) -> String {
        self.debug_info().unwrap_or_else(|_| "RustHLL()".to_string())
    }

    pub fn jaccard(&self, other: &RustHLL) -> PyResult<f64> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot compare sketches with different precisions"));
        }
        
        let mut intersection_registers = vec![0u8; self.registers.len()];
        let mut union_registers = vec![0u8; self.registers.len()];
        
        for i in 0..self.registers.len() {
            intersection_registers[i] = cmp::min(self.registers[i], other.registers[i]);
            union_registers[i] = cmp::max(self.registers[i], other.registers[i]);
        }
        
        let mut intersection_sketch = RustHLL::new(self.precision, None, None, None, None, None, None)?;
        intersection_sketch.registers = intersection_registers;
        intersection_sketch.alpha_mm = self.alpha_mm;
        
        let mut union_sketch = RustHLL::new(self.precision, None, None, None, None, None, None)?;
        union_sketch.registers = union_registers;
        union_sketch.alpha_mm = self.alpha_mm;
        
        let intersection = intersection_sketch.fast_mle_estimate();
        let union = union_sketch.fast_mle_estimate();
        
        if union == 0.0 {
            return Ok(0.0);
        }
        
        Ok(intersection / union)
    }

    pub fn is_empty(&self) -> PyResult<bool> {
        Ok(self.registers.iter().all(|&r| r == 0))
    }

    pub fn write(&self, filepath: &str) -> PyResult<()> {
        use std::fs::File;
        use std::io::Write;
        use bincode::serialize;

        let mut file = File::create(filepath)
            .map_err(|e| PyValueError::new_err(format!("Failed to create file: {}", e)))?;

        // Create a serializable struct with all the data
        #[derive(Serialize)]
        struct HLLData {
            registers: Vec<u8>,
            precision: u8,
            hash_size: u8,
            alpha_mm: f64,
            use_threading: bool,
            min_thread_batch: usize,
            compacted: bool,
        }

        let data = HLLData {
            registers: self.registers.clone(),
            precision: self.precision,
            hash_size: self.hash_size,
            alpha_mm: self.alpha_mm,
            use_threading: self.use_threading,
            min_thread_batch: self.min_thread_batch,
            compacted: self.compacted,
        };

        let serialized = serialize(&data)
            .map_err(|e| PyValueError::new_err(format!("Failed to serialize data: {}", e)))?;

        file.write_all(&serialized)
            .map_err(|e| PyValueError::new_err(format!("Failed to write to file: {}", e)))?;

        Ok(())
    }

    #[staticmethod]
    pub fn load(filepath: &str) -> PyResult<Self> {
        use std::fs::File;
        use std::io::Read;
        use bincode::deserialize;

        let mut file = File::open(filepath)
            .map_err(|e| PyValueError::new_err(format!("Failed to open file: {}", e)))?;

        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)
            .map_err(|e| PyValueError::new_err(format!("Failed to read file: {}", e)))?;

        #[derive(Deserialize)]
        struct HLLData {
            registers: Vec<u8>,
            precision: u8,
            hash_size: u8,
            alpha_mm: f64,
            use_threading: bool,
            min_thread_batch: usize,
            compacted: bool,
        }

        let data: HLLData = deserialize(&buffer)
            .map_err(|e| PyValueError::new_err(format!("Failed to deserialize data: {}", e)))?;

        Ok(RustHLL {
            registers: data.registers,
            precision: data.precision,
            mask: (1 << data.precision) - 1,
            alpha_mm: data.alpha_mm,
            use_threading: data.use_threading,
            min_thread_batch: data.min_thread_batch,
            compacted: data.compacted,
            hash_size: data.hash_size,
            hasher: RandomState::new(),
            memory_limit: None,
            debug: false,
            seed: None,
        })
    }

    pub fn set_memory_limit(&mut self, limit: Option<usize>) -> PyResult<()> {
        if let Some(limit) = limit {
            let current_memory = self.registers.len() * std::mem::size_of::<u8>();
            if current_memory > limit {
                return Err(PyValueError::new_err(format!(
                    "Current memory usage ({:.2} GB) exceeds new limit ({:.2} GB)",
                    current_memory as f64 / (1024.0 * 1024.0 * 1024.0),
                    limit as f64 / (1024.0 * 1024.0 * 1024.0)
                )));
            }
        }
        self.memory_limit = limit;
        Ok(())
    }

    pub fn get_memory_usage(&self) -> PyResult<usize> {
        Ok(self.registers.len() * std::mem::size_of::<u8>())
    }

    pub fn jaccard_register_match(&self, other: &RustHLL) -> PyResult<f64> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot compare sketches with different precisions"));
        }
        
        // Count matching and active registers
        let mut matching_count = 0;
        let mut active_count = 0;
        
        for i in 0..self.registers.len() {
            let val1 = self.registers[i];
            let val2 = other.registers[i];
            
            // Only consider registers where at least one sketch has a non-zero value
            if val1 > 0 || val2 > 0 {
                active_count += 1;
                if val1 == val2 {
                    matching_count += 1;
                }
            }
        }
        
        // If there are no active registers, return 0
        if active_count == 0 {
            return Ok(0.0);
        }
        
        // Jaccard is the proportion of matching registers among active ones
        Ok(matching_count as f64 / active_count as f64)
    }

    pub fn jaccard_minmax(&self, other: &RustHLL) -> PyResult<f64> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot compare sketches with different precisions"));
        }
        
        // Get cardinality estimates for both sketches
        let card_a = self.fast_mle_estimate();
        let card_b = other.fast_mle_estimate();
        
        // If either cardinality is 0, return 0 for Jaccard
        if card_a == 0.0 || card_b == 0.0 {
            return Ok(0.0);
        }
        
        // First get the register-based Jaccard estimate
        let register_jaccard = self.jaccard_register_match(other)?;
        
        // Calculate a weighted version using cardinality estimates
        // This helps correct for biases in the register-based approach
        let min_card = if card_a < card_b { card_a } else { card_b };
        let max_card = if card_a > card_b { card_a } else { card_b };
        
        // Scale the register Jaccard by cardinality ratio
        // This gives better estimates when sets have very different sizes
        let card_ratio = min_card / max_card;
        let weighted_jaccard = register_jaccard * (0.5 + 0.5 * card_ratio);
        
        // Ensure result is in [0,1]
        Ok(if weighted_jaccard > 1.0 { 1.0 } else { weighted_jaccard })
    }
    
    // Allow user to specify which Jaccard method to use
    pub fn jaccard_with_method(&self, other: &RustHLL, method: &str) -> PyResult<f64> {
        match method {
            "standard" => self.jaccard(other),
            "register_match" => self.jaccard_register_match(other),
            "minmax" => self.jaccard_minmax(other),
            "improved" => self.jaccard_improved(other),
            "minhash" => {
                // MinHash-like approach uses register encoding as hash value
                let mut matching_min_values = 0;
                let total_registers = self.registers.len() as f64;
                
                for i in 0..self.registers.len() {
                    if self.registers[i] > 0 && self.registers[i] == other.registers[i] {
                        matching_min_values += 1;
                    }
                }
                
                // Adjust for empty registers
                let empty_ratio = 1.0 - (matching_min_values as f64 / total_registers);
                Ok(1.0 - empty_ratio)
            },
            _ => Err(PyValueError::new_err(format!("Unknown Jaccard method: {}", method)))
        }
    }

    pub fn jaccard_improved(&self, other: &RustHLL) -> PyResult<f64> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot compare sketches with different precisions"));
        }
        
        // Get cardinality estimates
        let card_a = self.fast_mle_estimate();
        let card_b = other.fast_mle_estimate();
        
        // If either set is empty, Jaccard is 0
        if card_a == 0.0 || card_b == 0.0 {
            return Ok(0.0);
        }
        
        // First approach: Min/max registers (original)
        let mut intersection_registers = vec![0u8; self.registers.len()];
        let mut union_registers = vec![0u8; self.registers.len()];
        
        for i in 0..self.registers.len() {
            intersection_registers[i] = cmp::min(self.registers[i], other.registers[i]);
            union_registers[i] = cmp::max(self.registers[i], other.registers[i]);
        }
        
        let mut intersection_sketch = RustHLL::new(self.precision, None, None, None, None, None, None)?;
        intersection_sketch.registers = intersection_registers;
        intersection_sketch.alpha_mm = self.alpha_mm;
        
        let mut union_sketch = RustHLL::new(self.precision, None, None, None, None, None, None)?;
        union_sketch.registers = union_registers;
        union_sketch.alpha_mm = self.alpha_mm;
        
        let intersection_card = intersection_sketch.fast_mle_estimate();
        let union_card = union_sketch.fast_mle_estimate();
        
        // Second approach: register matching
        let mut matching_count = 0;
        let mut active_count = 0;
        
        for i in 0..self.registers.len() {
            let val1 = self.registers[i];
            let val2 = other.registers[i];
            
            if val1 > 0 || val2 > 0 {
                active_count += 1;
                if val1 == val2 {
                    matching_count += 1;
                }
            }
        }
        
        // Ensure we don't divide by zero
        let register_jaccard = if active_count > 0 {
            matching_count as f64 / active_count as f64
        } else {
            0.0
        };
        
        // Calculate set size difference
        let min_card = if card_a < card_b { card_a } else { card_b };
        let max_card = if card_a > card_b { card_a } else { card_b };
        let card_ratio = min_card / max_card;
        
        // Calculate the standard HLL Jaccard
        let jaccard_hll = if union_card > 0.0 { intersection_card / union_card } else { 0.0 };
        
        // Use different strategies based on the size ratio
        let jaccard = if card_ratio > 0.8 {
            // For similar size sets, traditional Jaccard works better
            jaccard_hll 
        } else if card_ratio > 0.5 {
            // For moderately different sizes, a mix of methods works best
            0.6 * jaccard_hll + 0.4 * register_jaccard
        } else {
            // For very different sized sets, register method is more reliable
            0.3 * jaccard_hll + 0.7 * register_jaccard
        };
        
        // Ensure result is in [0,1]
        Ok(jaccard)
    }

    #[getter]
    pub fn debug(&self) -> bool {
        self.debug
    }

    #[setter]
    pub fn set_debug(&mut self, value: bool) {
        self.debug = value;
    }

    #[getter]
    pub fn seed(&self) -> Option<u64> {
        self.seed
    }
}

fn get_index_and_count(hash: u64, precision: u8, hash_size: u8) -> (usize, u8) {
    // Get the index bits (lowest precision bits)
    let index = (hash & ((1 << precision) - 1)) as usize;
    
    // Get the remaining bits and count leading zeros
    let remaining = hash >> precision;
    let value = if remaining == 0 {
        hash_size - precision + 1
    } else {
        remaining.trailing_zeros() as u8 + 1
    };
    
    (index, value)
}

#[pymodule]
fn rust_hll(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<RustHLL>()?;
    Ok(())
} 