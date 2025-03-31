import numpy as np # type: ignore
import xxhash # type: ignore
from hammock.lib.abstractsketch import AbstractSketch
from typing import Dict, Optional, List

class HyperLogLog(AbstractSketch):
    def __init__(self, 
                 precision: int = 8, 
                 kmer_size: int = 0, 
                 window_size: int = 0, 
                 seed: Optional[int] = None,
                 debug: bool = False,
                 expected_cardinality: Optional[int] = None,
                 hash_size: int = 32):
        """Initialize HyperLogLog sketch.
        
        Args:
            precision: Number of bits for register indexing (4-{hash_size-1})
                      Lower values (8-10) work better for sparse sets
                      Higher values (14-24) work better for dense sets
            kmer_size: Size of k-mers (0 for whole string mode)
            window_size: Size of sliding window (0 or == kmer_size for no windowing)
            seed: Random seed for hashing
            debug: Whether to print debug information
            expected_cardinality: Expected number of unique items. If provided, precision will be adjusted.
            hash_size: Size of hash in bits (32 or 64)
        """
        super().__init__()
        
        if hash_size not in [32, 64]:
            raise ValueError("hash_size must be 32 or 64")
        
        if expected_cardinality is not None:
            # Adjust precision based on expected cardinality
            if expected_cardinality < 1000:
                precision = 8  # For small sets, use lower precision
            elif expected_cardinality < 10000:
                precision = 12  # For medium sets
            elif expected_cardinality < 100000:
                precision = 18  # For larger sets
            else:
                precision = 22  # For very large sets
        
        if precision < 4:
            raise ValueError(f"Precision must be at least 4")
        if precision >= hash_size:
            raise ValueError(f"Precision must be less than hash_size ({hash_size})")
        
        if kmer_size < 0:
            raise ValueError("k-mer size must be non-negative")
        
        if window_size and window_size < kmer_size:
            raise ValueError("Window size must be >= kmer size")
        
        self.precision = precision
        self.num_registers = 1 << precision
        self.registers = np.zeros(self.num_registers, dtype=np.int8)
        self.kmer_size = kmer_size
        self.window_size = window_size if window_size else kmer_size
        self.debug = debug
        self.seed = seed if seed is not None else 42
        self.hash_size = hash_size

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
    
    @staticmethod
    def _hash32_int(x: int, seed: int = 0) -> int:
        """32-bit hash function for integers.
        
        Args:
            x: Integer value to hash
            seed: Random seed for hashing
            
        Returns:
            32-bit hash value as integer
        """
        hasher = xxhash.xxh32(seed=seed)
        hasher.update(x.to_bytes(8, byteorder='little'))
        return hasher.intdigest()

    # def hash64(self, x: int) -> int:
    #     """Instance method to hash an integer using the instance's seed."""
    #     return self._hash64_int(x, seed=self.seed)

    @staticmethod
    def _hash_str(s: bytes, seed: int = 0, hash_size: int = 32) -> int:
        """Hash a string using xxhash.
        
        Args:
            s: String to hash
            seed: Random seed for hashing
            hash_size: Size of hash in bits (32 or 64)
            
        Returns:
            Hash value as integer
        """
        if hash_size == 32:
            hasher = xxhash.xxh32(seed=seed)
        else:
            hasher = xxhash.xxh64(seed=seed)
        hasher.update(s)
        return hasher.intdigest()

    def hash_str(self, s: bytes) -> int:
        """Instance method to hash a string using the instance's seed."""
        return self._hash_str(s, seed=self.seed, hash_size=self.hash_size)

    def _rho(self, hash_val: int) -> int:
        """Calculate position of leftmost 1-bit."""
        hash_val >>= self.precision
        pos = 1
        while (hash_val & 1) == 0 and pos <= self.hash_size - self.precision:
            pos += 1
            hash_val >>= 1
        return pos

    def _process_kmer(self, kmer: str) -> None:
        """Process a single k-mer or string."""
        hash_val = self.hash_str(kmer.encode())
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)
    
    def _register_counts(self) -> np.ndarray:
        """Get counts of registers."""
        return np.bincount(self.registers, minlength=self.hash_size - self.precision + 1)

    def add_string_with_minimizers(self, s: str) -> None:
        """Add a string to the sketch using minimizer windowing scheme.
        
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
                h = self._hash_str(kmer.encode(), seed=self.seed, hash_size=self.hash_size)
                if h < min_hash:
                    min_hash = h
                    min_pos = j
            
            # Process the minimizer
            self._process_kmer(window[min_pos:min_pos + self.kmer_size])

    def add_string(self, s: str) -> None:
        """Add a string to the sketch.
        
        Processes all k-mers in the string.
        If kmer_size is 0, processes whole string.
        
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

        # Process all k-mers
        for i in range(len(s) - self.kmer_size + 1):
            self._process_kmer(s[i:i + self.kmer_size])

    def add_batch(self, strings: List[str]) -> None:
        """Add multiple strings to the sketch.
        
        Args:
            strings: List of strings to add to the sketch
        """
        for s in strings:
            self.add_string(s)

    def hash64_int(self, x: int) -> int:
        """Instance method to hash an integer using the instance's seed."""
        return self._hash64_int(x, seed=self.seed)
    
    def hash32_int(self, x: int) -> int:
        """Instance method to hash an integer using the instance's seed."""
        return self._hash32_int(x, seed=self.seed)

    def add_int(self, value: int) -> None:
        """Add an integer to the sketch."""
        self.item_count += 1
        if self.hash_size == 32:
            hash_val = self._hash32_int(value, self.seed)
        else:
            hash_val = self._hash64_int(value, self.seed)
        idx = hash_val & (self.num_registers - 1)
        rank = self._rho(hash_val)
        self.registers[idx] = max(self.registers[idx], rank)

    def estimate_cardinality(self, method: str = 'ertl_mle') -> float:
        """Estimate the cardinality of the multiset."""
        self.register_counts = self._register_counts()
        if np.all(self.registers == 0):
            return 0.0
        if method == 'original':
            return self.original_estimate()
        elif method == 'ertl_improved':
            return self.ertl_improved_estimate()
        elif method == 'ertl_mle':
            return self.ertl_mle_estimate()
        else:
            raise ValueError(f"Invalid method: {method}")
        
    def original_estimate(self) -> float:
        m = float(self.num_registers)
        alpha = self._get_alpha(m)
        
        register_harmonics = np.power(2.0, -self.registers)
        raw_estimate = alpha * m * m / np.sum(register_harmonics)
        
        # Small range correction
        if raw_estimate <= 2.5 * m:
            v = np.sum(self.registers == 0)
            if v > 0:
                raw_estimate = m * np.log(m / float(v))
                return raw_estimate
        
        # Large range correction
        max_value = pow(2, self.hash_size - self.precision)
        if raw_estimate > max_value / 32.0:
            # Calculate the argument for logarithm with safety check
            log_arg = 1.0 - raw_estimate / max_value
            # If the argument is close to or less than zero, cap the estimate
            if log_arg <= 0.000001:  # Small epsilon to avoid numerical issues
                return 0.9 * max_value  # Return slightly less than max_value
            
            # Safe to calculate the logarithm now
            raw_estimate = -max_value * np.log(log_arg)
        
        return raw_estimate

    def _get_tau(self, x) -> float:
        """Calculate tau value for HyperLogLog bias correction.
        Args:
            x: Input value between 0 and 1
        Returns:
            Bias correction factor tau(x)
        """
        if x == 0.0 or x == 1.0:
            return 0.0
        
        # Initialize variables
        z = 1.0 - x  # z starts as complement of x
        tmp = 0.0
        y = 1.0
        zprev = x       # zprev (z previous) starts as x
        
        # Iterate until convergence (z stops changing)
        while zprev != z:
            x = x ** 0.5    # Square root of x
            zprev = z          # Save previous z
            y *= 0.5       # Halve y each iteration (geometric scaling)
            tmp = 1.0 - x   # Calculate temporary value
            z -= tmp * tmp * y  # Update z
            
        return z / 3.0
            
    def ertl_improved_estimate(self) -> float:
        """Estimate cardinality using Ertl's improved method."""
        norm_const = 0.7213  # 1/(2*ln(2))
        m = float(self.num_registers)
        # counts = self.register_counts()
    
        # Get counts of maxed registers
        n_maxed_registers = self.register_counts[self.hash_size - self.precision + 1]
        non_maxreg_frac = (m - n_maxed_registers)/m
         # Get counts of zero registers
        zero_reg_frac = self.register_counts[0]/m
    
        # Apply bias corrections
        tau = self._get_tau(non_maxreg_frac)
        z = m * tau
        # Process intermediate registers
        for i in range(self.hash_size - self.precision, 0, -1):
            z += self.register_counts[i]  # Add count for this register value
            z *= 0.5 # geometric scaling
        
        #bias correction for zero registers
        sigma = self._get_sigma(zero_reg_frac)
        z+= m*sigma
        #bias correction for non-maxed registers
        return m * norm_const * m / z

    
    def _get_probs(self, est: float) -> float:
        """Calculate probabilities for each register value.
        
        Args:
            est: Estimated cardinality
        Returns:
            Total probability
        """
        total_prob = 0.0
        
        # Find optimization bounds
        # k_min: Find first non-zero register count to skip empty low registers
        # This improves performance and numerical stability since very low registers
        # often have no counts
        k_min = 0
        while k_min < 64 and self.register_counts[k_min] == 0:
            k_min += 1
        k_min = max(1, k_min)  # Ensure we start at least at 1 to avoid numerical issues with 0
        
        # k_max: Find last non-zero register count to skip empty high registers
        # This avoids unnecessary computation for register values that can't occur
        # given our precision parameter
        k_max = self.hash_size - self.precision
        while k_max > 0 and self.register_counts[k_max] == 0:
            k_max -= 1
        k_max = min(self.hash_size - self.precision, k_max)  # Ensure we don't exceed maximum possible value
        
        # Only calculate probabilities for non-zero register counts within our bounds
        # This optimization:
        # 1. Reduces computation by skipping empty registers
        # 2. Improves numerical stability by avoiding extreme probability values
        # 3. Focuses calculation on the registers that actually contribute to the estimate
        for j in range(k_min, k_max + 1):
            if self.register_counts[j] > 0:
                # Calculate probability of getting j leading zeros
                # p = P(value = j) = P(zeros >= j) - P(zeros >= j+1)
                p = (1.0 - 2.0 ** -j) ** est - (1.0 - 2.0 ** -(j+1)) ** est
                total_prob += p * self.register_counts[j]
            
        return total_prob

    def ertl_mle_estimate(self, relative_error: float = 1e-2) -> float:
        """Estimate cardinality using Ertl's Maximum Likelihood Estimation method."""
        num_registers = 1 << self.precision
        max_register_value = self.hash_size - self.precision
        
        # Check if all registers are maxed out - fix potential index out of bounds
        register_counts_size = len(self.register_counts)
        if max_register_value + 1 < register_counts_size and self.register_counts[max_register_value + 1] == num_registers:
            return float('inf')
        
        # Find bounds for non-zero register values
        min_nonzero_value = 0
        while min_nonzero_value < self.hash_size and self.register_counts[min_nonzero_value] == 0:
            min_nonzero_value += 1
        min_value_bound = max(1, min_nonzero_value)  # Ensure minimum of 1 for numerical stability
        
        max_nonzero_value = min(max_register_value, register_counts_size - 1)
        while max_nonzero_value > 0 and self.register_counts[max_nonzero_value] == 0:
            max_nonzero_value -= 1
        max_value_bound = min(max_register_value, max_nonzero_value)
        
        # Calculate initial harmonic mean with proper scaling
        harmonic_mean = 0.0
        for register_value in range(max_value_bound, min_value_bound - 1, -1):
            harmonic_mean = 0.5 * harmonic_mean + self.register_counts[register_value]
        harmonic_mean = harmonic_mean * (2.0 ** -min_value_bound)  # Scale by power of 2
        
        # Initialize estimation parameters
        maxed_register_count = 0
        if max_register_value + 1 < register_counts_size:
            maxed_register_count = self.register_counts[max_register_value + 1]
        if max_register_value < register_counts_size and max_value_bound < register_counts_size:
            maxed_register_count += self.register_counts[max_value_bound]
        
        zero_correction = harmonic_mean + self.register_counts[0]
        active_registers = num_registers - self.register_counts[0]
        
        # Initial estimate using improved bound
        prev_estimate = harmonic_mean
        if max_register_value + 1 < register_counts_size:
            prev_estimate += self.register_counts[max_register_value + 1] * (2.0 ** -max_register_value)
        
        if prev_estimate <= 1.5 * zero_correction:
            cardinality_estimate = active_registers / (0.5 * prev_estimate + zero_correction)
        else:
            cardinality_estimate = (active_registers / prev_estimate) * np.log1p(prev_estimate / zero_correction)
        
        # Adjust relative error based on number of registers
        relative_error /= np.sqrt(num_registers)
        estimate_change = cardinality_estimate
        
        # Main estimation loop with improved convergence
        while estimate_change > cardinality_estimate * relative_error:
            # Calculate binary exponent for scaling
            binary_exponent = int(np.floor(np.log2(cardinality_estimate))) + 1
            scaled_estimate = cardinality_estimate * (2.0 ** -max(max_value_bound + 1, binary_exponent + 2))
            scaled_estimate_squared = scaled_estimate * scaled_estimate
            
            # Taylor series approximation for probability function
            prob_estimate = (scaled_estimate - 
                            scaled_estimate_squared/3 + 
                            (scaled_estimate_squared * scaled_estimate_squared) * 
                            (1/45 - scaled_estimate_squared/472.5))
            
            # Update probability estimate for each register value
            total_probability = maxed_register_count * prob_estimate
            for register_value in range(max_value_bound - 1, min_value_bound - 1, -1):
                if register_value < register_counts_size:  # Check bounds
                    prob_complement = 1.0 - prob_estimate
                    prob_estimate = ((scaled_estimate + prob_estimate * prob_complement) / 
                                   (scaled_estimate + prob_complement))
                    scaled_estimate *= 2
                    if self.register_counts[register_value] > 0:
                        total_probability += self.register_counts[register_value] * prob_estimate
                else:
                    # Skip this iteration if register_value is out of bounds
                    prob_complement = 1.0 - prob_estimate
                    prob_estimate = ((scaled_estimate + prob_estimate * prob_complement) / 
                                   (scaled_estimate + prob_complement))
                    scaled_estimate *= 2
            
            total_probability += cardinality_estimate * zero_correction
            
            # Update estimate using secant method
            if prev_estimate < total_probability <= active_registers:
                estimate_change *= ((total_probability - active_registers) / 
                                  (prev_estimate - total_probability))
            else:
                estimate_change = 0
            
            prev_estimate += estimate_change
        
        return prev_estimate * num_registers

 
    def _get_sigma(self, x: float) -> float:
        """Calculate sigma value for HyperLogLog bias correction.
        
        Args:
            x: Input value between 0 and 1 (fraction of empty registers)
            
        Returns:
            Bias correction factor sigma(x)
        """
        if x == 1.0:
            return float('inf')
        
        z = x
        zprev = 0.0  # z previous
        y = 1.0
        
        while z != zprev:
            x = x * x    
            zprev = z 
            z += x * y  # Update z
            y += y      # Double y
            
            if np.isnan(z):
                # If we hit a numerical instability, return last valid value
                return zprev
        return z

    def _get_alpha(self, m: float) -> float:
        """Get alpha constant based on number of registers."""
        if m == 16:
            return 0.673
        elif m == 32:
            return 0.697
        elif m == 64:
            return 0.709
        else:
            # Add a very small correction factor to reduce overestimation
            return 0.7213 / (1.0 + 1.079 / m) * 0.985

    def raw_estimate(self) -> float:
        """Calculate the raw cardinality estimate before corrections.
        
        Returns:
            Raw estimated cardinality
        """
        # Sum 2^(-max(register)) for each register
        # sum_inv = sum(math.pow(2.0, -max(0, x)) for x in self.registers)
        registers = np.array(self.registers, dtype=(np.float32 if self.hash_size == 32 else np.float64))
        sum_inv = np.sum(np.exp2(-registers))
        
        # Calculate alpha_m * m^2 / sum
        alpha_m = self.get_alpha()
        estimate = alpha_m * (self.num_registers * self.num_registers) / sum_inv
        
        return estimate

    def merge(self, other: 'HyperLogLog') -> None:
        """Merge another HLL sketch into this one.
        
        This modifies the current sketch by taking the element-wise maximum of 
        its registers
        with the other sketch's registers, effectively combining their 
        cardinality estimates.
        
        Args:
            other: Another HyperLogLog sketch to merge into this one
            
        Raises:
            ValueError: If the sketches have different precision values
        """
        if not isinstance(other, HyperLogLog):
            raise TypeError("Can only merge with another HyperLogLog sketch")
        if self.precision != other.precision:
            raise ValueError("Cannot merge HyperLogLog sketches with different precisions")
        if self.kmer_size != other.kmer_size:
            raise ValueError("Cannot merge HyperLogLog sketches with different k-mer sizes")
        
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
        if not isinstance(other, HyperLogLog):
            raise TypeError("Can only compare with another HyperLogLog sketch")
        if self.precision != other.precision:
            raise ValueError("Cannot compare HLLs with different precision")
        if self.seed != other.seed:
            raise ValueError("HyperLogLogs must have same seed for comparison")
        
        if self.is_empty() or other.is_empty():
            return 0.0        
        return self.estimate_jaccard_registers(other)  # Use register method by default

    def estimate_jaccard_registers(self, other: 'HyperLogLog') -> float:
        """Estimate Jaccard similarity with another HyperLogLog using register min/max."""
        
        # Return 0 if either sketch is empty

        
        # Get non-zero masks for both register sets
        nonzero1 = self.registers != 0
        nonzero2 = other.registers != 0
        
        # Only consider registers where at least one sketch has a non-zero value
        active_registers = nonzero1 | nonzero2
        if not active_registers.any():
            return 0.0
        
        # Count matching non-zero registers
        matching_nonzero = (self.registers[active_registers] == other.registers
        [active_registers]).sum()
        total_active = active_registers.sum()
        
        # Jaccard is the proportion of matching registers among active ones
        return matching_nonzero / total_active
        
        # Create union and intersection sketches with same seed
        #union = HyperLogLog(precision=self.precision, seed=self.seed)
        intersection = HyperLogLog(precision=self.precision, seed=self.seed)
        # intersection = self.estimate_intersection(other)
        # union = self.estimate_union(other)
        # # Compute register-wise max (union) and min (intersection)
        # np.maximum(self.registers, other.registers, out=union.registers)
        # np.minimum(self.registers, other.registers, out=intersection.registers)
        # # for i in range(self.num_registers):
        # #     union.registers[i] = max(self.registers[i], other.registers[i])
        # #     intersection.registers[i] = min(self.registers[i], other.registers
        # #     [i])
        
        # # Estimate Jaccard similarity
        # union_card = union.estimate_cardinality()
        # if union_card == 0:
        #     return 0.0
        
        # intersection_card = intersection.estimate_cardinality()
        # return intersection_card / union_card

    def estimate_jaccard_iep(self, other: 'HyperLogLog') -> float:
        """Estimate Jaccard similarity using inclusion-exclusion principle."""
        # Return 0 if either sketch is empty
        if self.is_empty() or other.is_empty():
            return 0.0
        
        # Get cardinality estimates
        card_a = self.estimate_cardinality()
        card_b = other.estimate_cardinality()
        
        # Create union sketch
        card_union = self.estimate_union(other)
        
        # Use inclusion-exclusion principle
        card_intersection = max(0.0, card_a + card_b - card_union)
        
        if card_union == 0:
            return 0.0
        
        return card_intersection / card_union

    def estimate_jaccard_original(self, other: 'HyperLogLog') -> float:
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
            seed=np.array([self.seed]),
            hash_size=np.array([self.hash_size])
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
            seed=int(data['seed'][0]),
            hash_size=int(data['hash_size'][0])
        )
        sketch.registers = data['registers']
        return sketch

    def similarity_values(self, other: 'AbstractSketch') -> Dict[str, float]:
        """Calculate similarity values using HyperLogLog.
        
        Returns:
            Dictionary containing 'jaccard_similarity'
        """
        if not isinstance(other, HyperLogLog):
            raise ValueError("Can only compare with another HyperLogLog sketch")
        if self.kmer_size != other.kmer_size:
            raise ValueError(f"Cannot compare HyperLogLog sketches with different k-mer sizes ({self.kmer_size} vs {other.kmer_size})")
        
        jaccard = self.estimate_jaccard(other)
        return {'jaccard_similarity': jaccard}

    def is_empty(self) -> bool:
        """Check if sketch is empty."""
        return np.all(self.registers == 0)
       