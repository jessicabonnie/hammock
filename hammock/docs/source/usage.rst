Usage Guide
===========

Basic Usage
----------

The Hammock library provides three main functions for comparing files:

1. ``compare_bed_files()`` for BED file comparison
2. ``compare_sequence_files()`` for sequence file comparison
3. ``compare_files()`` for general file comparison

Comparing BED Files
-----------------

.. code-block:: python

    from hammock import compare_bed_files

    # Basic usage
    results = compare_bed_files("file1.bed", "file2.bed")
    print(results)

    # Advanced usage with custom parameters
    results = compare_bed_files(
        "file1.bed",
        "file2.bed",
        sketch_type="hyperloglog",
        precision=12,
        subA=0.8,
        subB=1.0,
        expA=0.5,
        use_rust=True
    )

Comparing Sequence Files
----------------------

.. code-block:: python

    from hammock import compare_sequence_files

    # Basic usage
    results = compare_sequence_files("file1.fa", "file2.fa")
    print(results)

    # Advanced usage with custom parameters
    results = compare_sequence_files(
        "file1.fa",
        "file2.fa",
        kmer_size=10,
        window_size=50,
        sketch_type="minhash",
        num_hashes=128
    )

General File Comparison
---------------------

.. code-block:: python

    from hammock import compare_files

    # Compare any two files
    results = compare_files(
        "file1.txt",
        "file2.txt",
        mode='A',  # Compare intervals only
        sketch_type="hyperloglog",
        precision=12,
        num_hashes=64,
        subA=1.0,
        subB=1.0,
        expA=0.5,
        use_rust=True
    )

Command Line Interface
--------------------

Hammock also provides a command-line interface:

.. code-block:: bash

    # Compare BED files
    hammock compare-bed file1.bed file2.bed

    # Compare sequence files
    hammock compare-seq file1.fa file2.fa

    # Compare any files with custom parameters
    hammock compare file1.txt file2.txt --mode A --sketch-type hyperloglog --precision 12

Parameters
---------

Common Parameters
~~~~~~~~~~~~~~~

* ``sketch_type``: Type of sketch to use
    * "hyperloglog" (default)
    * "minhash"
    * "minimizer"
    * "exact"
* ``precision``: Precision for HyperLogLog sketching (default: 12)
* ``num_hashes``: Number of hashes for MinHash sketching (default: 64)

BED File Parameters
~~~~~~~~~~~~~~~~~

* ``subA``: Subsampling rate for intervals (0 to 1, default: 1.0)
* ``subB``: Subsampling rate for points (0 to 1, default: 1.0)
* ``expA``: Power of 10 exponent for A-type intervals (default: 0.5)
* ``use_rust``: Whether to use Rust implementation for HyperLogLog (default: True)

Sequence File Parameters
~~~~~~~~~~~~~~~~~~~~~~

* ``kmer_size``: Size of k-mers for sequence sketching (default: 8)
* ``window_size``: Size of sliding window for sequence sketching (default: 40)

Output Format
------------

The comparison functions return a dictionary containing similarity metrics:

.. code-block:: python

    {
        'jaccard': 0.75,  # Jaccard similarity
        'containment': 0.8,  # Containment similarity
        'intersection': 1000,  # Size of intersection
        'union': 2000  # Size of union
    }

Best Practices
-------------

1. **Memory Management**
   * For large files, consider using subsampling (``subA`` and ``subB`` parameters)
   * Use the Rust implementation for better performance with HyperLogLog

2. **Accuracy vs. Speed**
   * Higher precision values give more accurate results but use more memory
   * More hash functions in MinHash give more accurate results but take longer

3. **File Formats**
   * BED files should be properly formatted with at least 3 columns (chr, start, end)
   * Sequence files should be in FASTA format

4. **Performance Optimization**
   * Use the appropriate mode for your data type
   * Consider using the Rust implementation for large datasets
   * Adjust k-mer and window sizes based on your sequence data 