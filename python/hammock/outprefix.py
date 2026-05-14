from typing import Optional


def get_new_prefix(outprefix: str,
                   sketch_type: str,
                   mode: str,
                   num_hashes: Optional[int] = None,
                   precision: Optional[int] = None,
                   subA: Optional[float] = None,
                   subB: Optional[float] = None,
                   expA: Optional[float] = None,
                   kmer_size: Optional[int] = None,
                   window_size: Optional[int] = None) -> str:
    if sketch_type == "minhash":
        outprefix = f"{outprefix}_mh_n{num_hashes}_jacc{mode}"
    elif sketch_type == "hyperloglog":
        outprefix = f"{outprefix}_hll_p{precision}_jacc{mode}"
    elif sketch_type == "minimizer":
        outprefix = f"{outprefix}_mnmzr_p{precision}_jacc{mode}"
    else:
        outprefix = f"{outprefix}_exact_jacc{mode}"

    if expA is not None and expA > 0:
        outprefix = f"{outprefix}_expA{expA:.2f}"
    elif subA is not None and subA != 1.0:
        outprefix = f"{outprefix}_A{subA:.2f}"

    if subB is not None and subB != 1.0:
        # 4 decimals when subB < 0.01 so we don't collide on rounding (e.g.
        # 0.0001 vs 0.001 both -> "0.00" under :.2f). Keep :.2f at the
        # historical resolution for subB >= 0.01 so existing filenames stay
        # stable.
        fmt = ".4f" if subB < 0.01 else ".2f"
        outprefix = f"{outprefix}_B{format(subB, fmt)}"

    if mode == "D" and kmer_size is not None and window_size is not None:
        outprefix = f"{outprefix}_k{kmer_size}_w{window_size}"

    return outprefix
