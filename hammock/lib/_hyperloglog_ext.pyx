# cython_hll_batch.pyx
import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint64_t, uint32_t, uint8_t, int8_t
from libc.string cimport memcmp

# Import xxhash for compatibility
import xxhash

# Type definitions
ctypedef np.int8_t DTYPE_t
ctypedef np.uint64_t HASH_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef class CythonHLLBatch:
    """Cython-optimized HyperLogLog batch operations."""
    
    cdef:
        int precision
        int num_registers
        int kmer_size
        int window_size
        int hash_size
        uint64_t seed
        np.ndarray registers
        int8_t[:] registers_view
        
    def __init__(self, int precision, int kmer_size=0, int window_size=0, 
                 uint64_t seed=42, int hash_size=32):
        self.precision = precision
        self.num_registers = 1 << precision
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size > 0 else kmer_size
        self.hash_size = hash_size
        self.seed = seed
        
        # Initialize registers array with memory view for fast access
        self.registers = np.zeros(self.num_registers, dtype=np.int8)
        self.registers_view = self.registers
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef uint64_t _hash_bytes(self, bytes data):
        """Fast hash function using xxhash."""
        if self.hash_size == 32:
            return xxhash.xxh32(data, seed=self.seed).intdigest()
        else:
            return xxhash.xxh64(data, seed=self.seed).intdigest()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int _rho(self, uint64_t hash_val):
        """Calculate position of leftmost 1-bit (leading zeros + 1)."""
        hash_val >>= self.precision
        cdef int pos = 1
        while (hash_val & 1) == 0 and pos <= self.hash_size - self.precision:
            pos += 1
            hash_val >>= 1
        return pos
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _process_kmer_fast(self, bytes kmer):
        """Fast k-mer processing with minimal Python overhead."""
        cdef uint64_t hash_val = self._hash_bytes(kmer)
        cdef int idx = hash_val & (self.num_registers - 1)
        cdef int rank = self._rho(hash_val)
        
        # Update register with maximum value
        if rank > self.registers_view[idx]:
            self.registers_view[idx] = rank
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_string_batch(self, list strings):
        """Batch process multiple strings with optimized loops."""
        cdef:
            int i, j, k
            int str_len, num_kmers
            bytes string_bytes, kmer_bytes
            str string
            
        # Process each string in the batch
        for i in range(len(strings)):
            string = strings[i]
            string_bytes = string.encode('utf-8')
            str_len = len(string_bytes)
            
            if self.kmer_size == 0:
                # Whole string mode
                self._process_kmer_fast(string_bytes)
            else:
                # K-mer mode
                if str_len >= self.kmer_size:
                    num_kmers = str_len - self.kmer_size + 1
                    for j in range(num_kmers):
                        kmer_bytes = string_bytes[j:j + self.kmer_size]
                        self._process_kmer_fast(kmer_bytes)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_string_batch_with_minimizers(self, list strings):
        """Batch process with minimizer windowing scheme."""
        cdef:
            int i, j, k
            int str_len, num_windows
            bytes string_bytes, window_bytes, kmer_bytes
            str string
            uint64_t min_hash, current_hash
            int min_pos
            
        for i in range(len(strings)):
            string = strings[i]
            string_bytes = string.encode('utf-8')
            str_len = len(string_bytes)
            
            if str_len < self.kmer_size:
                continue
                
            if self.kmer_size == 0:
                # Whole string mode
                self._process_kmer_fast(string_bytes)
            elif self.window_size == self.kmer_size:
                # K-mer mode without windowing
                num_windows = str_len - self.kmer_size + 1
                for j in range(num_windows):
                    kmer_bytes = string_bytes[j:j + self.kmer_size]
                    self._process_kmer_fast(kmer_bytes)
            else:
                # Windowing mode with minimizers
                num_windows = str_len - self.window_size + 1
                for j in range(num_windows):
                    window_bytes = string_bytes[j:j + self.window_size]
                    min_hash = <uint64_t>-1  # Max uint64_t
                    min_pos = 0
                    
                    # Find minimum hash in this window
                    for k in range(self.window_size - self.kmer_size + 1):
                        kmer_bytes = window_bytes[k:k + self.kmer_size]
                        current_hash = self._hash_bytes(kmer_bytes)
                        if current_hash < min_hash:
                            min_hash = current_hash
                            min_pos = k
                    
                    # Process the minimizer
                    kmer_bytes = window_bytes[min_pos:min_pos + self.kmer_size]
                    self._process_kmer_fast(kmer_bytes)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def add_kmers_batch(self, list kmers):
        """Directly process a batch of k-mers (pre-extracted)."""
        cdef:
            int i
            bytes kmer_bytes
            str kmer
            
        for i in range(len(kmers)):
            kmer = kmers[i]
            kmer_bytes = kmer.encode('utf-8')
            self._process_kmer_fast(kmer_bytes)
    
    def get_registers(self):
        """Return the current register values."""
        return np.asarray(self.registers_view)
    
    def merge_registers(self, np.ndarray[DTYPE_t, ndim=1] other_registers):
        """Merge another set of registers into this one."""
        cdef:
            int i
            int8_t[:] other_view = other_registers
            
        for i in range(self.num_registers):
            if other_view[i] > self.registers_view[i]:
                self.registers_view[i] = other_view[i] 