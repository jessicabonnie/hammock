Welcome to Hammock's documentation!
================================

Hammock is a Python library for efficiently comparing files, particularly BED files and sequence files, using various sketching algorithms. It provides a simple and intuitive API for computing similarity metrics between files while maintaining memory efficiency.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   api
   examples
   contributing

Features
--------

* Efficient comparison of BED files and sequence files
* Multiple sketching algorithms (HyperLogLog, MinHash, Minimizer)
* Rust implementation for improved performance
* Memory-efficient processing of large files
* Simple and intuitive API
* Comprehensive documentation and examples

Quick Start
----------

Install Hammock:

.. code-block:: bash

    pip install hammock

Compare two BED files:

.. code-block:: python

    import hammock

    # Compare two BED files
    results = hammock.compare_bed_files('file1.bed', 'file2.bed')
    print(f"Jaccard similarity: {results['jaccard']}")
    print(f"Containment: {results['containment']}")

Compare two sequence files:

.. code-block:: python

    # Compare two sequence files
    results = hammock.compare_sequence_files('seq1.fa', 'seq2.fa')
    print(f"Sequence similarity: {results['jaccard']}")

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 