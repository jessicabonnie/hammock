"""Reference-genome resolution for the bed2fasta feature (Mode D BED→FASTA).

**v1 policy: no network access during a `hammock` run.** `resolve_reference`
only ever returns a path to an already-local, indexed FASTA:

  * a local FASTA path (``.fa``/``.fasta``/``.fna``, optionally ``.gz``) is used
    directly — this is the "bring your own reference" case;
  * a keyword (``hg38``, ``mm10``, …) is looked up in ``cache_dir``. If it isn't
    cached yet we raise with the exact ``hammock fetch-ref`` command to populate
    it (run once on a networked login node — HPC compute nodes are firewalled).

All downloading lives in `fetch_reference`, reached only via the explicit
``hammock fetch-ref`` subcommand — never during a comparison run.

The cache layout for a keyword ``K`` is::

    <cache_dir>/K.fa        decompressed FASTA
    <cache_dir>/K.fa.fai    faidx index
    <cache_dir>/K.fa.done   sentinel written last (guards partial fetches)
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
from urllib.parse import urlparse

# Built-in keyword → UCSC goldenPath bigZips ``.fa.gz`` URL map. Soft-masked
# assemblies are fine: digest's minimizer hashing is case-insensitive.
_KEYWORD_URLS = {
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "mm10": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
    "mm39": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
    "hs1":  "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz",
}

# Friendly aliases → canonical keyword.
_KEYWORD_ALIASES = {
    "grch38": "hg38", "gencode38": "hg38",
    "grch37": "hg19",
    "chm13": "hs1", "t2t": "hs1", "chm13v2": "hs1", "chm13v2.0": "hs1",
}

_FASTA_EXTS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")

# Refuse to decompress anything larger than this (guards zip-bombs / runaway
# fetches from a bad URL). T2T/large mammalian genomes fit comfortably.
_MAX_DECOMPRESSED_BYTES = 20 * 1024 * 1024 * 1024  # 20 GiB


def default_cache_dir() -> str:
    """Where cached references live unless ``--ref-cache-dir`` overrides it."""
    return os.environ.get(
        "HAMMOCK_REF_CACHE",
        os.path.join(os.path.expanduser("~"), ".hammock", "refs"),
    )


def canonical_keyword(spec: str) -> str | None:
    """Return the canonical keyword for ``spec`` (via alias table), or None."""
    key = spec.strip().lower()
    key = _KEYWORD_ALIASES.get(key, key)
    return key if key in _KEYWORD_URLS else None


def known_keywords() -> list[str]:
    return sorted(set(_KEYWORD_URLS) | set(_KEYWORD_ALIASES))


def _looks_like_url(spec: str) -> bool:
    return spec.startswith(("http://", "https://"))


def _is_gz(path: str) -> bool:
    return path.lower().endswith(".gz")


def _cache_fasta_path(cache_dir: str, keyword: str) -> str:
    return os.path.join(cache_dir, f"{keyword}.fa")


def _is_complete(fasta: str) -> bool:
    """A cached reference is usable iff its ``.done`` sentinel and a fresh
    ``.fai`` both exist (index no older than the FASTA)."""
    done = fasta + ".done"
    fai = fasta + ".fai"
    if not (os.path.exists(done) and os.path.exists(fai) and os.path.exists(fasta)):
        return False
    try:
        return os.path.getmtime(fai) >= os.path.getmtime(fasta)
    except OSError:
        return False


# --------------------------------------------------------------------------
# Resolution (run path — no network)
# --------------------------------------------------------------------------

def resolve_reference(spec: str, cache_dir: str | None = None) -> str:
    """Resolve ``spec`` to a local, indexed FASTA path. Never downloads.

    Raises ValueError with actionable guidance if a keyword isn't cached or the
    spec is unrecognised.
    """
    if not spec:
        raise ValueError("Empty reference specification.")
    cache_dir = cache_dir or default_cache_dir()

    # 1. Local FASTA path — bring your own reference.
    if os.path.exists(spec):
        low = spec.lower()
        if not low.endswith(_FASTA_EXTS):
            raise ValueError(
                f"Reference path '{spec}' does not look like a FASTA "
                f"(expected one of {_FASTA_EXTS})."
            )
        if _is_gz(spec):
            # Decompress+index into the cache once (bedtools getfasta wants a
            # plain, faidx-able FASTA). Keyed by the source path.
            key = _safe_key("local_" + os.path.abspath(spec))
            target = _cache_fasta_path(cache_dir, key)
            if _is_complete(target):
                return target
            os.makedirs(cache_dir, exist_ok=True)
            _publish_from_gz_file(spec, target)
            return target
        return _ensure_indexed_local(spec, cache_dir)

    # 2. Keyword → cache lookup (no download here).
    keyword = canonical_keyword(spec)
    if keyword is not None:
        target = _cache_fasta_path(cache_dir, keyword)
        if _is_complete(target):
            return target
        raise ValueError(
            f"Reference keyword '{spec}' (→ {keyword}) is not cached in "
            f"{cache_dir}. Fetch it once on a networked (login) node:\n"
            f"    hammock fetch-ref {keyword} --ref-cache-dir {cache_dir}\n"
            f"(HPC compute nodes are typically firewalled, so hammock never "
            f"downloads during a run.)"
        )

    # 3. A URL: resolve against the cache only (never download in a run).
    if _looks_like_url(spec):
        target = _cache_fasta_path(cache_dir, _safe_key("url_" + spec))
        if _is_complete(target):
            return target
        raise ValueError(
            f"URL reference '{spec}' is not cached in {cache_dir}. Fetch it once "
            f"on a networked (login) node:\n"
            f"    hammock fetch-ref {spec} --ref-cache-dir {cache_dir}\n"
            f"then re-run with the same --ref/--ref-cache-dir."
        )

    raise ValueError(
        f"Unrecognised reference '{spec}'. Provide a local FASTA path, or a "
        f"known keyword: {', '.join(known_keywords())}."
    )


def _ensure_indexed_local(fasta: str, cache_dir: str) -> str:
    """Ensure ``<fasta>.fai`` exists. If the FASTA's directory is read-only
    (the normal shared ``/data`` case), symlink it into ``cache_dir`` and index
    there. Returns the path bedtools should use as ``-fi``."""
    if _fai_fresh(fasta):
        return fasta

    fasta_dir = os.path.dirname(os.path.abspath(fasta)) or "."
    if os.access(fasta_dir, os.W_OK):
        # Best-effort: build now with samtools if present, else let bedtools
        # getfasta create the .fai at extraction time (dir is writable).
        _faidx(fasta, required=False)
        return fasta

    # Read-only source dir → index a symlink in the cache instead.
    os.makedirs(cache_dir, exist_ok=True)
    link = os.path.join(cache_dir, _safe_key("local_" + os.path.abspath(fasta))
                        + os.path.splitext(fasta)[1])
    if not os.path.islink(link):
        # os.symlink can't atomically replace; remove a stale non-symlink.
        if os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(fasta), link)
    if not (os.path.exists(link + ".fai") and _fai_fresh(link)):
        _faidx(link)
    return link


def _fai_fresh(fasta: str) -> bool:
    fai = fasta + ".fai"
    try:
        return os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(fasta)
    except OSError:
        return False


# --------------------------------------------------------------------------
# Fetch path (explicit `hammock fetch-ref` — network allowed)
# --------------------------------------------------------------------------

def fetch_reference(spec: str, cache_dir: str | None = None,
                    force: bool = False) -> str:
    """Download + decompress + index a reference into ``cache_dir``. Returns the
    local FASTA path. Only invoked by the ``hammock fetch-ref`` subcommand."""
    cache_dir = cache_dir or default_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)

    keyword = canonical_keyword(spec)
    if keyword is not None:
        url = _KEYWORD_URLS[keyword]
        target = _cache_fasta_path(cache_dir, keyword)
    elif _looks_like_url(spec):
        scheme = urlparse(spec).scheme
        if scheme not in ("http", "https"):
            raise ValueError(f"Only http(s) URLs are allowed; got scheme '{scheme}'.")
        url = spec
        target = _cache_fasta_path(cache_dir, _safe_key("url_" + spec))
    elif os.path.exists(spec):
        # Already local: just index it (or decompress a .gz into the cache).
        return resolve_reference(spec, cache_dir)
    else:
        raise ValueError(
            f"'{spec}' is neither a known keyword, an http(s) URL, nor a local "
            f"path. Known keywords: {', '.join(known_keywords())}."
        )

    if _is_complete(target) and not force:
        print(f"[hammock] already cached: {target}", file=sys.stderr)
        return target

    print(f"[hammock] fetching {url}\n           → {target}", file=sys.stderr)
    _publish_from_url(url, target)
    print(f"[hammock] done: {target}", file=sys.stderr)
    return target


def _publish_from_url(url: str, target: str) -> None:
    """Download ``url`` (a ``.fa.gz``) and atomically publish an indexed FASTA
    at ``target``. Never leaves a partial file at the canonical path."""
    import tempfile
    import urllib.request

    tmp_gz = None
    try:
        fd, tmp_gz = tempfile.mkstemp(dir=os.path.dirname(target) or ".", suffix=".gz.tmp")
        os.close(fd)
        req = urllib.request.Request(url, headers={"User-Agent": "hammock-fetch-ref"})
        with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_gz, "wb") as out:
            shutil.copyfileobj(resp, out, length=1024 * 1024)
        with open(tmp_gz, "rb") as fh:
            if fh.read(2) != b"\x1f\x8b":
                raise ValueError(
                    f"Downloaded file from {url} is not gzip (redirect to an "
                    f"error page?). Aborting.")
        _publish_from_gz_file(tmp_gz, target, source_url=url)
    finally:
        if tmp_gz and os.path.exists(tmp_gz):
            os.remove(tmp_gz)


def _publish_from_gz_file(gz_path: str, target: str, source_url: str = "") -> None:
    """Decompress ``gz_path`` → ``target`` (atomic), index, write ``.done``."""
    import tempfile

    tmp_fa = None
    try:
        fd, tmp_fa = tempfile.mkstemp(dir=os.path.dirname(target) or ".", suffix=".fa.tmp")
        written = 0
        with gzip.open(gz_path, "rb") as src, os.fdopen(fd, "wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                written += len(chunk)
                if written > _MAX_DECOMPRESSED_BYTES:
                    raise ValueError(
                        f"Decompressed size exceeds {_MAX_DECOMPRESSED_BYTES} "
                        f"bytes; refusing (possible zip-bomb).")
                dst.write(chunk)
        with open(tmp_fa, "rb") as fh:
            if not fh.read(1).startswith(b">"):
                raise ValueError("Decompressed file does not look like FASTA "
                                 "(no leading '>').")
        os.replace(tmp_fa, target)
        tmp_fa = None
        _faidx(target)
        _write_done(target, source_url, os.path.getsize(target))
    finally:
        if tmp_fa and os.path.exists(tmp_fa):
            os.remove(tmp_fa)


def _write_done(fasta: str, source_url: str, size: int) -> None:
    import datetime
    import tempfile
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(fasta) or ".", suffix=".done.tmp")
    with os.fdopen(fd, "w") as fh:
        fh.write(f"source={source_url}\nsize={size}\n"
                 f"date={datetime.datetime.now().isoformat()}\n")
    os.replace(tmp, fasta + ".done")


def _faidx(fasta: str, required: bool = True) -> None:
    """Build ``<fasta>.fai`` with samtools.

    With ``required=True`` (the publish path, which writes a ``.done`` sentinel
    afterwards) samtools is mandatory — otherwise ``.done`` would claim a
    complete cache with no index. With ``required=False`` (a writable
    bring-your-own FASTA) a missing samtools is tolerated: bedtools getfasta
    will build the ``.fai`` itself at extraction time."""
    samtools = shutil.which("samtools")
    if samtools is None:
        if required:
            raise RuntimeError(
                "samtools not found, so the reference cannot be indexed. Load it "
                "('ml samtools') and retry.")
        return  # bedtools getfasta will create the .fai itself.
    subprocess.run([samtools, "faidx", fasta], check=True,
                   capture_output=True, text=True)


def _safe_key(s: str) -> str:
    """A filesystem-safe cache key derived from an arbitrary string (no path
    traversal). Used for URLs and local-path-derived cache entries."""
    import hashlib
    return "cache_" + hashlib.sha256(s.encode()).hexdigest()[:16]
