import numpy as np # type: ignore
import xxhash # type: ignore
from hammock.lib.abstractsketch import AbstractSketch

class HyperLogLog(AbstractSketch):
    def __init__(self, 
                 precision: int = 8, 
                 kmer_size: int = 0, 
                 window_size: int = 0, 
                 seed: int = 0,
                 debug: bool = False):
        """Initialize HyperLogLog sketch.
        
        Args:
            precision: Number of bits for register indexing (4-16)
                      Lower values (8-10) work better for sparse sets
                      Higher values (14-16) work better for dense sets
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
        """
        if not 4 <= precision <= 16:
            raise ValueError("Precision must be between 4 and 16")
        
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.seed = seed
        self.registers = np.zeros(self.num_registers, dtype=np.uint8)
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.debug = debug

        self.item_count = 0
        
        # Calculate alpha_mm (bias correction factor)
        if self.num_registers == 16:
            self.alpha_mm = 0.673
        elif self.num_registers == 32:
            self.alpha_mm = 0.697
        elif self.num_registers == 64:
            self.alpha_mm = 0.709
        else:
            self.alpha_mm = 0.7213 / (1 + 1.079 / self.num_registers)


    @staticmethod
    def _hash64_int(x: int, seed: int = 0) -> int:
        """64-bit hash function for integers.
        
        Args:
            x: Integer value to hash
            seed: Random seed for hashing
            
        Returns:
            64-bit hash value as integer
        """
        hasher = xxhash.xxh64(seed=seed)
        hasher.update(x.to_bytes(8, byteorder='little'))
        return hasher.intdigest()

    # def hash64(self, x: int) -> int:
    #     """Instance method to hash an integer using the instance's seed."""
    #     return self._hash64_int(x, seed=self.seed)

    @staticmethod
    def _hash_str(s: bytes, seed: int = 0) -> int:
        """Hash a string using xxhash."""
        hasher = xxhash.xxh64(seed=seed)
        hasher.update(s)
        return hasher.intdigest()

    def hash_str(self, s: bytes) -> int:
        """Instance method to hash a string using the instance's seed."""
        return self._hash_str(s, seed=self.seed)

    def _rho(self, hash_val: int) -> int:
        """Calculate position of leftmost 1-bit."""
        hash_val >>= self.precision
        pos = 1
        while (hash_val & 1) == 0 and pos <= 64:
            pos += 1
            hash_val >>= 1
        return pos

    def _process_kmer(self, kmer: str) -> None:
        """Process a single k-mer or string."""
        hash_val = self._hash_str(kmer.encode())
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)

    def add_string(self, s: str) -> None:
        """Add a string to the sketch.
        
        Processes the string according to kmer_size and window_size settings.
        If kmer_size is 0, processes whole string.
        If window_size equals kmer_size, processes all k-mers.
        Otherwise uses minimizer scheme within windows.
        
        Args:
            s: String to add to sketch
        """
        if len(s) < self.kmer_size:
            return
        
        self.item_count += 1

        # Whole string mode
        if self.kmer_size == 0:
            self._process_kmer(s)
            return

        # K-mer mode without windowing
        if self.window_size == self.kmer_size:
            for i in range(len(s) - self.kmer_size + 1):
                self._process_kmer(s[i:i + self.kmer_size])
            return

        # Windowing mode
        for i in range(len(s) - self.window_size + 1):
            window = s[i:i + self.window_size]
            min_hash = float('inf')
            min_pos = 0
            
            # Find minimum hash in this window
            for j in range(self.window_size - self.kmer_size + 1):
                kmer = window[j:j + self.kmer_size]
                h = self._hash_str(kmer.encode())
                if h < min_hash:
                    min_hash = h
                    min_pos = j
            
            # Process the minimizer
            self._process_kmer(window[min_pos:min_pos + self.kmer_size])

    def hash64_int(self, x: int) -> int:
        """Instance method to hash an integer using the instance's seed."""
        return self._hash64_int(x, seed=self.seed)

    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        self.item_count += 1
        hash_val = self.hash64_int(value)
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)

    def estimate_cardinality(self) -> float:
        """Estimate the cardinality of the set using HyperLogLog algorithm."""
        registers = np.array(self.registers, dtype=np.float64)
        sum_inv = np.sum(np.exp2(-registers))
        estimate = self.alpha_mm * (self.num_registers ** 2) / sum_inv
        if self.debug:
            print(f"DEBUG cardinality: estimate={estimate:.1f}, items={self.item_count}")
        return estimate

    def raw_estimate(self) -> float:
        """Calculate the raw cardinality estimate before corrections.
        
        Returns:
            Raw estimated cardinality
        """
        # Sum 2^(-max(register)) for each register
        # sum_inv = sum(math.pow(2.0, -max(0, x)) for x in self.registers)
        registers = np.array(self.registers, dtype=np.float64)
        sum_inv = np.sum(np.exp2(-registers))
        
        # Calculate alpha_m * m^2 / sum
        alpha_m = self.get_alpha()
        estimate = alpha_m * (self.num_registers * self.num_registers) / sum_inv
        
        return estimate

    def merge(self, other: 'HyperLogLog') -> None:
        """Merge another HLL sketch into this one.
        
        This modifies the current sketch by taking the element-wise maximum of its registers
        with the other sketch's registers, effectively combining their cardinality estimates.
        
        Args:
            other: Another HyperLogLog sketch to merge into this one
            
        Raises:
            ValueError: If the sketches have different precision values
        """
        if self.precision != other.precision:
            raise ValueError("Cannot merge HLLs with different precision")
        
        # Take element-wise maximum and modify self.registers in-place
        np.maximum(self.registers, other.registers, out=self.registers)

    def estimate_intersection(self, other: 'HyperLogLog') -> float:
        """Estimate intersection cardinality with another HLL.
        
        Uses inclusion-exclusion principle with union estimate.
        
        Args:
            other: Another HyperLogLog sketch
            
        Returns:
            Estimated size of intersection between sketches
            
        Raises:
            ValueError: If sketches have different precision values
        """

        if self.precision != other.precision:
            raise ValueError("Cannot compute intersection of HLLs with different precision")
        
        # Basic inclusion-exclusion
        a = self.estimate_cardinality()
        b = other.estimate_cardinality()
        union = self.estimate_union(other)
        intersection = max(0.0, a + b - union)
        
        if self.debug:
            print(f"DEBUG: items={self.item_count}/{other.item_count}, a={a:.1f}, b={b:.1f}, union={union:.1f}, intersection={intersection:.1f}")
        return intersection

    def estimate_jaccard(self, other: 'HyperLogLog') -> float:
        """Estimate Jaccard similarity with another HyperLogLog sketch."""
        if self.precision != other.precision:
            raise ValueError("Cannot compare HLLs with different precision")
        
        intersection = self.estimate_intersection(other)
        union = self.estimate_union(other)
        
        if union == 0:
            return 0.0
        
        jaccard = intersection / union
        if self.debug:
            print(f"DEBUG: jaccard={jaccard:.3f}")
        return jaccard

    def estimate_union(self, other: 'HyperLogLog') -> float:
        """Estimate union cardinality with another HyperLogLog sketch."""
        if self.precision != other.precision:
            raise ValueError("Cannot compute union of HLLs with different precision")
        
        # Simply take max of registers and estimate
        merged = HyperLogLog(self.precision, self.kmer_size, self.window_size, self.seed, self.debug)
        merged.registers = np.maximum(self.registers, other.registers)
        return merged.estimate_cardinality()

    def get_alpha(self) -> float:
        """Get alpha correction factor based on number of registers.
        
        The alpha factor corrects for bias in the HyperLogLog algorithm.
        Values are based on the original HyperLogLog paper.
        
        Returns:
            Alpha correction factor as a float
        """
        # This is redundant since alpha_mm is already calculated in __init__
        return self.alpha_mm

    def write(self, filepath: str) -> None:
        """Write sketch to file in binary format.
        
        Args:
            filepath: Path to output file
        """
        np.savez_compressed(
            filepath,
            registers=self.registers,
            precision=np.array([self.precision]),
            kmer_size=np.array([self.kmer_size]),
            window_size=np.array([self.window_size]),
            seed=np.array([self.seed])
        )

    @classmethod
    def load(cls, filepath: str) -> 'HyperLogLog':
        """Load sketch from file in binary format.
        
        Args:
            filepath: Path to input file
            
        Returns:
            HyperLogLog object loaded from file
        """
        data = np.load(filepath)
        sketch = cls(
            precision=int(data['precision'][0]),
            kmer_size=int(data['kmer_size'][0]),
            window_size=int(data['window_size'][0]),
            seed=int(data['seed'][0])
        )
        sketch.registers = data['registers']
        return sketch