use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::cmp;

/// A fast HyperLogLog implementation in Rust
#[pyclass]
struct RustHLL {
    registers: Vec<u8>,
    precision: usize,
    mask: u64,
    alpha_mm: f64,
}

#[pymethods]
impl RustHLL {
    /// Create a new HyperLogLog sketch with the given precision
    #[new]
    fn new(precision: usize) -> PyResult<Self> {
        if !(4..=16).contains(&precision) {
            return Err(PyValueError::new_err("Precision must be between 4 and 16"));
        }
        
        // Number of registers (m = 2^precision)
        let m = 1 << precision;
        
        // Pre-compute constants
        let alpha = match m {
            16 => 0.673,
            32 => 0.697,
            64 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / m as f64),
        };
        
        Ok(RustHLL {
            registers: vec![0; m],
            precision,
            mask: (m as u64) - 1,
            alpha_mm: alpha * (m as f64) * (m as f64),
        })
    }
    
    /// Add a string value to the sketch
    fn add(&mut self, value: &str) -> PyResult<()> {
        let (index, count) = self.get_index_and_count(value);
        self.registers[index] = cmp::max(self.registers[index], count);
        Ok(())
    }
    
    /// Add a batch of string values to the sketch
    fn add_batch(&mut self, values: Vec<&str>) -> PyResult<()> {
        for value in values {
            self.add(value)?;
        }
        Ok(())
    }
    
    /// Estimate the cardinality (number of unique elements)
    fn estimate(&self) -> f64 {
        // Compute the harmonic mean
        let mut sum = 0.0;
        for &value in &self.registers {
            sum += 2.0f64.powi(-(value as i32));
        }
        
        let estimate = self.alpha_mm / sum;
        let m = self.registers.len() as f64;
        
        // Apply correction for small cardinalities
        if estimate <= 2.5 * m {
            // Count number of registers equal to 0
            let v = self.registers.iter().filter(|&&x| x == 0).count() as f64;
            if v > 0.0 {
                return m * (m / v).ln();
            }
        }
        
        // Apply correction for large cardinalities
        if estimate > u32::MAX as f64 / 30.0 {
            return -(u32::MAX as f64) * (1.0 - estimate / (u32::MAX as f64)).ln();
        }
        
        estimate
    }
    
    /// Merge another HyperLogLog sketch into this one
    fn merge(&mut self, other: &RustHLL) -> PyResult<()> {
        if self.precision != other.precision {
            return Err(PyValueError::new_err("Cannot merge sketches with different precisions"));
        }
        
        for i in 0..self.registers.len() {
            self.registers[i] = cmp::max(self.registers[i], other.registers[i]);
        }
        
        Ok(())
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
        
        // Create a merged sketch
        let mut merged = RustHLL::new(self.precision)?;
        merged.merge(self)?;
        merged.merge(other)?;
        
        // Estimate the Jaccard similarity using the inclusion-exclusion principle
        let self_est = self.estimate();
        let other_est = other.estimate();
        let union_est = merged.estimate();
        
        if union_est == 0.0 {
            return Ok(1.0); // Both are empty, consider them identical
        }
        
        // Jaccard = |A ∩ B| / |A ∪ B| = (|A| + |B| - |A ∪ B|) / |A ∪ B|
        let intersection_est = self_est + other_est - union_est;
        let similarity = intersection_est / union_est;
        
        // Clamp to valid range [0, 1]
        Ok(similarity.max(0.0).min(1.0))
    }
}

impl RustHLL {
    /// Hash the value and return the index and bit pattern
    fn get_index_and_count(&self, value: &str) -> (usize, u8) {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        let hash = hasher.finish();
        
        // Use the first precision bits as the register index
        let index = (hash & self.mask) as usize;
        
        // Use the position of the leftmost 1-bit in the remaining bits for the count
        let bits = hash >> self.precision;
        let leading_zeros = bits.leading_zeros() as u8 + 1;
        
        (index, leading_zeros)
    }
}

#[pymodule]
fn rust_hll(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<RustHLL>()?;
    Ok(())
} 