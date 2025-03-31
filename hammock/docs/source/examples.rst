Examples
========

Basic Examples
-------------

Comparing BED Files
~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_bed_files

    # Compare two BED files
    results = compare_bed_files("peaks1.bed", "peaks2.bed")
    print(f"Jaccard similarity: {results['jaccard']}")
    print(f"Containment: {results['containment']}")

    # Compare with custom parameters
    results = compare_bed_files(
        "peaks1.bed",
        "peaks2.bed",
        sketch_type="hyperloglog",
        precision=14,  # Higher precision for more accuracy
        subA=0.8,  # Subsample intervals
        use_rust=True  # Use Rust implementation
    )

Comparing Sequence Files
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_sequence_files

    # Compare two sequence files
    results = compare_sequence_files("seq1.fa", "seq2.fa")
    print(f"Jaccard similarity: {results['jaccard']}")

    # Compare with custom parameters
    results = compare_sequence_files(
        "seq1.fa",
        "seq2.fa",
        kmer_size=10,  # Larger k-mers
        window_size=50,  # Larger windows
        sketch_type="minhash",
        num_hashes=128  # More hashes for better accuracy
    )

Advanced Examples
---------------

Batch Processing
~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_bed_files
    import glob
    import pandas as pd

    # Compare multiple BED files
    bed_files = glob.glob("peaks/*.bed")
    results = []

    for i, file1 in enumerate(bed_files):
        for file2 in bed_files[i+1:]:
            result = compare_bed_files(file1, file2)
            results.append({
                'file1': file1,
                'file2': file2,
                'jaccard': result['jaccard'],
                'containment': result['containment']
            })

    # Convert to DataFrame
    df = pd.DataFrame(results)
    print(df)

Memory-Efficient Processing
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_bed_files

    # Process large files with subsampling
    results = compare_bed_files(
        "large_peaks1.bed",
        "large_peaks2.bed",
        subA=0.5,  # Subsample 50% of intervals
        subB=0.5,  # Subsample 50% of points
        sketch_type="hyperloglog",
        precision=12
    )

Custom Comparison Modes
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_files

    # Compare both intervals and points
    results = compare_files(
        "file1.txt",
        "file2.txt",
        mode='C',  # Compare both intervals and points
        sketch_type="hyperloglog",
        precision=12
    )

    # Compare points only
    results = compare_files(
        "file1.txt",
        "file2.txt",
        mode='B',  # Compare points only
        sketch_type="minhash",
        num_hashes=64
    )

Real-World Examples
-----------------

ChIP-seq Peak Comparison
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_bed_files
    import pandas as pd

    # Compare ChIP-seq peaks from different conditions
    results = compare_bed_files(
        "control_peaks.bed",
        "treatment_peaks.bed",
        sketch_type="hyperloglog",
        precision=14,  # High precision for accurate comparison
        subA=0.8,  # Subsample to handle large peak sets
        use_rust=True  # Use Rust for better performance
    )

    print(f"Peak overlap similarity: {results['jaccard']}")
    print(f"Treatment peaks contained in control: {results['containment']}")

Sequence Similarity Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from hammock import compare_sequence_files
    import glob

    # Compare multiple sequence files
    fasta_files = glob.glob("sequences/*.fa")
    results = []

    for i, file1 in enumerate(fasta_files):
        for file2 in fasta_files[i+1:]:
            result = compare_sequence_files(
                file1,
                file2,
                kmer_size=10,  # Larger k-mers for more specific matching
                window_size=50,  # Larger windows for better context
                sketch_type="minhash",
                num_hashes=128  # More hashes for better accuracy
            )
            results.append({
                'file1': file1,
                'file2': file2,
                'similarity': result['jaccard']
            })

    # Print results
    for r in results:
        print(f"{r['file1']} vs {r['file2']}: {r['similarity']:.3f}") 