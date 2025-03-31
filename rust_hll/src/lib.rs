use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
#[allow(unused_imports)]
use std::collections::hash_map::DefaultHasher;
#[allow(unused_imports)]
use std::hash::{Hash, Hasher};
use std::cmp;
use std::thread;
use std::sync::{Arc, Mutex};
use num_cpus;
use bincode;
use serde::{Serialize, Deserialize};

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
}

// Internal implementation
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
                let (idx, val) = get_index_and_count(hashes[i], self.precision, self.hash_size);
                indices[i] = idx;
                values[i] = val;
            }
            
            // Second pass: update registers in sorted order for better cache locality
            let mut idx_val_pairs: Vec<(usize, u8)> = indices.into_iter()
                .zip(values.into_iter())
                .collect();
            
            // Sort by index for better cache locality
            idx_val_pairs.sort_by_key(|&(idx, _)| idx);
            
            // Update registers in sorted order
            for (idx, val) in idx_val_pairs {
                self.update_register(idx, val);
            }
        } else {
            // For small batches, process directly
            for &hash in hashes {
                let (idx, val) = get_index_and_count(hash, self.precision, self.hash_size);
                self.update_register(idx, val);
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
            if max_val > 30 {
                // Compact registers if max value is too large
                for r in &mut self.registers {
                    if *r > 30 {
                        *r = 30;
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
    #[pyo3(signature = (precision, use_threading=None, min_thread_batch=None, hash_size=None))]
    pub fn new(precision: u8, use_threading: Option<bool>, min_thread_batch: Option<usize>, hash_size: Option<u8>) -> PyResult<Self> {
        // Validate hash_size
        let hash_size = hash_size.unwrap_or(32);  // Default to 32-bit hashing
        if hash_size != 32 && hash_size != 64 {
            return Err(PyValueError::new_err("hash_size must be 32 or 64"));
        }

        // Check if precision would exceed the hash length
        if precision >= hash_size {
            return Err(PyValueError::new_err(format!("Precision must be < {} ({}-bit hash)", hash_size, hash_size)));
        }

        let m = 1 << precision;
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
            min_thread_batch: min_thread_batch.unwrap_or(10000),
            compacted: false,
            hash_size,
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
        let mut hasher = DefaultHasher::new();
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
            let mut hasher = DefaultHasher::new();
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
        let num_threads = num_cpus::get();
        let chunk_size = (values.len() + num_threads - 1) / num_threads;
        
        let values_arc = Arc::new(values);
        let registers_arc = Arc::new(Mutex::new(vec![0u8; self.registers.len()]));
        let precision = self.precision;
        let hash_size = self.hash_size;  // Clone hash_size before the loop
        
        let mut handles = vec![];
        
        for i in 0..num_threads {
            let start = i * chunk_size;
            let end = cmp::min(start + chunk_size, values_arc.len());
            
            if start >= end {
                break;
            }
            
            let values_clone = Arc::clone(&values_arc);
            let registers_clone = Arc::clone(&registers_arc);
            
            let handle = thread::spawn(move || {
                let mut local_registers = vec![0u8; 1 << precision];
                
                for value in &values_clone[start..end] {
                    let mut hasher = DefaultHasher::new();
                    value.hash(&mut hasher);
                    let hash = if hash_size == 32 {
                        hasher.finish() as u32 as u64
                    } else {
                        hasher.finish()
                    };
                    
                    let (idx, val) = get_index_and_count(hash, precision, hash_size);
                    if local_registers[idx] < val {
                        local_registers[idx] = val;
                    }
                }
                
                let mut registers = registers_clone.lock().unwrap();
                for i in 0..registers.len() {
                    if registers[i] < local_registers[i] {
                        registers[i] = local_registers[i];
                    }
                }
            });
            
            handles.push(handle);
        }
        
        for handle in handles {
            handle.join().unwrap();
        }
        
        let final_registers = registers_arc.lock().unwrap();
        self.registers.copy_from_slice(&final_registers);
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
        
        let mut intersection_sketch = RustHLL::new(self.precision, None, None, None)?;
        intersection_sketch.registers = intersection_registers;
        intersection_sketch.alpha_mm = self.alpha_mm;
        
        let mut union_sketch = RustHLL::new(self.precision, None, None, None)?;
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
        })
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