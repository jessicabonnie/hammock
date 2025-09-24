# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

"""
High-performance BED parser with fast ingestion into C++ HyperLogLog.

Exports a minimal API expected by IntervalSketch:
 - CppBedParser(...)
 - parse_file(path, debug=False) -> list[str]
 - parse_into_hll(path, fast_hll, debug=False, expA: int = 0) -> None

Notes:
 - Uses xxhash64 for stable 64-bit hashing when inserting via add_hash64.
 - Keeps Python I/O for portability; heavy work (subsampling, per-base iteration,
   hashing) is in Cython to minimize Python overhead.
"""

from libc.stdint cimport uint64_t
import xxhash

cdef inline uint64_t _hash64_bytes(bytes b, uint64_t seed):
    # xxhash Python object is not GIL-free; call under GIL.
    return <uint64_t>xxhash.xxh64(b, seed=seed).intdigest()

cdef inline bint _is_header_or_blank(str line) noexcept:
    if not line:
        return True
    if line[0] == '#':
        return True
    # quick header detection by first token
    try:
        first = line.split()[0].lower()
    except Exception:
        return True
    return first in ('chromosome', 'start', 'end', 'chrom', 'chromosome_name')

cdef class CppBedParser:
    cdef:
        int precision
        int hash_size
        uint64_t seed
        str mode
        double subsample_a
        double subsample_b
        str separator
    
    def __init__(self, int precision=12, int hash_size=64, uint64_t seed=42, 
                 mode='A', double subsample_a=1.0, double subsample_b=1.0,
                 separator='\t'):
        if mode not in ('A','B','C'):
            raise ValueError("CppBedParser: mode must be A, B, or C")
        if hash_size not in (32, 64):
            raise ValueError("hash_size must be 32 or 64")
        if subsample_a < 0.0 or subsample_a > 1.0 or subsample_b < 0.0 or subsample_b > 1.0:
            raise ValueError("subsample rates must be in [0,1]")
        self.precision = precision
        self.hash_size = hash_size
        self.seed = seed
        self.mode = mode
        self.subsample_a = subsample_a
        self.subsample_b = subsample_b
        self.separator = separator

    def parse_file(self, path: str, debug: bool=False):
        """Return strings representing intervals/points according to mode.
        Warning: For large BEDs this can be memory-heavy; prefer parse_into_hll.
        """
        cdef list out = []
        cdef str line
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if _is_header_or_blank(line):
                    continue
                cols = line.split()
                if len(cols) < 3:
                    continue
                chrval = cols[0]
                try:
                    start = int(cols[1])
                    end = int(cols[2])
                except Exception:
                    continue
                if start < 0 or end <= start:
                    continue

                if self.mode in ('A','C'):
                    if self.mode == 'C' and self.subsample_a < 1.0:
                        # subsample intervals deterministically via hash of interval string
                        istring = f"{chrval}{self.separator}{start}{self.separator}{end}{self.separator}A"
                        h = xxhash.xxh32(istring, seed=31337).intdigest()
                        thr = int(self.subsample_a * ((1<<32) - 1))
                        if h <= thr:
                            out.append(istring)
                    else:
                        out.append(f"{chrval}{self.separator}{start}{self.separator}{end}{self.separator}A")

                if self.mode in ('B','C'):
                    if self.subsample_b >= 1.0:
                        for x in range(start, end):
                            out.append(f"{chrval}{self.separator}{x}")
                    else:
                        thr = int(self.subsample_b * ((1<<32) - 1))
                        for x in range(start, end):
                            pstr = f"{chrval}{self.separator}{x}"
                            h32 = xxhash.xxh32(pstr, seed=31337).intdigest()
                            if h32 <= thr:
                                out.append(pstr)
        if debug:
            print(f"CppBedParser.parse_file produced {len(out)} strings")
        return out

    def parse_into_hll(self, path: str, fast_hll, debug: bool=False, int expA=0):
        """Stream parse BED file and insert hashed values directly into FastHyperLogLog.
        Uses xxhash64 to produce stable 64-bit hashes and calls fast_hll.add_hash64.
        """
        cdef str line
        cdef uint64_t h64
        cdef int x
        cdef int i
        cdef bint include_interval
        cdef bint do_A = self.mode in ('A','C')
        cdef bint do_B = self.mode in ('B','C')
        cdef int mult = 0
        if expA > 0:
            mult = 1
            for i in range(expA):
                mult *= 10
        cdef int base_interval_count = 0
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if _is_header_or_blank(line):
                    continue
                cols = line.split()
                if len(cols) < 3:
                    continue
                chrval = cols[0]
                try:
                    start = int(cols[1]); end = int(cols[2])
                except Exception:
                    continue
                if start < 0 or end <= start:
                    continue

                if do_A:
                    istring = f"{chrval}{self.separator}{start}{self.separator}{end}{self.separator}A".encode('utf-8')
                    include_interval = True
                    if self.mode == 'C' and self.subsample_a < 1.0:
                        h32 = xxhash.xxh32(istring, seed=31337).intdigest()
                        thrA = int(self.subsample_a * ((1<<32) - 1))
                        include_interval = (h32 <= thrA)

                    if include_interval:
                        base_interval_count += 1
                        if mult > 0:
                            # Add exactly mult total copies, each uniquely suffixed
                            for i in range(mult):
                                istring2 = istring + str(i+1).encode('utf-8')
                                h64 = _hash64_bytes(istring2, <uint64_t>self.seed)
                                fast_hll.add_hash64(h64)
                        else:
                            # No multiplicity: add the base interval once
                            h64 = _hash64_bytes(istring, <uint64_t>self.seed)
                            fast_hll.add_hash64(h64)

                if do_B:
                    # Subsampling for points is only valid in mode C
                    if not (self.mode == 'C' and self.subsample_b < 1.0):
                        for x in range(start, end):
                            pbytes = f"{chrval}{self.separator}{x}".encode('utf-8')
                            h64 = _hash64_bytes(pbytes, <uint64_t>self.seed)
                            fast_hll.add_hash64(h64)
                    else:
                        thr = int(self.subsample_b * ((1<<32) - 1))
                        for x in range(start, end):
                            pstr = f"{chrval}{self.separator}{x}"
                            h32 = xxhash.xxh32(pstr, seed=31337).intdigest()
                            if h32 <= thr:
                                pbytes = pstr.encode('utf-8')
                                h64 = _hash64_bytes(pbytes, <uint64_t>self.seed)
                                fast_hll.add_hash64(h64)
        if debug:
            print("CppBedParser.parse_into_hll completed; intervals added:", base_interval_count)
        return base_interval_count
 
