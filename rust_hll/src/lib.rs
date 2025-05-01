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
}

// Internal implementation
impl RustHLL {
    /// Process a batch of hashes in chunks for better cache locality and memory efficiency
    fn process_hash_batch(&mut self, hashes: &[u64]) {
        // Process in smaller chunks to limit memory usage
        const CHUNK_SIZE: usize = 256;  // Increased back to original value
        
        // Pre-allocate hash vector with fixed capacity
        let mut hash_buffer = Vec::with_capacity(CHUNK_SIZE);
            
        for chunk in hashes.chunks(CHUNK_SIZE) {
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
    #[pyo3(signature = (precision, use_threading=None, min_thread_batch=None, hash_size=None, memory_limit=None, debug=None))]
    pub fn new(precision: u8, use_threading: Option<bool>, min_thread_batch: Option<usize>, hash_size: Option<u8>, memory_limit: Option<usize>, debug: Option<bool>) -> PyResult<Self> {
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
            if debug.unwrap_or(false) {
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
            
            if debug.unwrap_or(false) {
                println!("Rust memory check:");
                println!("  - Estimated memory: {} bytes", estimated_memory);
                println!("  - Estimated with buffer: {} bytes", estimated_memory_with_buffer);
                println!("  - Memory limit: {} bytes", limit);
            }
            
            if estimated_memory_with_buffer > limit {
                return Err(PyValueError::new_err(format!(
                    "Estimated memory usage with buffer ({:.2} GB) exceeds limit ({:.2} GB)",
                    estimated_memory_with_buffer as f64 / (1024.0 * 1024.0 * 1024.0),
                    limit as f64 / (1024.0 * 1024.0 * 1024.0)
                )));
            }
        }

        let alpha_mm = match precision {
            4 => 0.673,
            5 => 0.697,
            6 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / m as f64),
        };

        Ok(RustHLL {
            registers: vec![0; m],
            precision,
            mask: (1 << precision) - 1,
            alpha_mm,
            use_threading: use_threading.unwrap_or(true),
            min_thread_batch: min_thread_batch.unwrap_or(1000),
            compacted: false,
            hash_size,
            hasher: RandomState::new(),
            memory_limit,
            debug: debug.unwrap_or(false),
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
        
        let mut intersection_sketch = RustHLL::new(self.precision, None, None, None, None, None)?;
        intersection_sketch.registers = intersection_registers;
        intersection_sketch.alpha_mm = self.alpha_mm;
        
        let mut union_sketch = RustHLL::new(self.precision, None, None, None, None, None)?;
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