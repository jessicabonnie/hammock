Installation
============

Basic Installation
-----------------

You can install Hammock using pip:

.. code-block:: bash

    pip install hammock

Development Installation
-----------------------

For development installation:

.. code-block:: bash

    git clone https://github.com/jbonnie1/hammock.git
    cd hammock
    pip install -e ".[dev]"

Dependencies
------------

Hammock requires the following Python packages:

* numpy >= 1.19.0
* pandas >= 1.2.0
* scipy >= 1.6.0

Optional Dependencies
--------------------

For development, the following additional packages are installed:

* pytest >= 6.0 (for testing)
* black >= 21.0 (for code formatting)
* mypy >= 0.800 (for type checking)
* flake8 >= 3.9 (for linting)
* sphinx >= 4.0.0 (for documentation)

Rust Dependencies
----------------

For optimal performance with HyperLogLog sketching, Hammock can use a Rust implementation. To use this feature:

1. Install Rust from https://www.rust-lang.org/tools/install
2. The Rust implementation will be automatically compiled during installation

System Requirements
------------------

* Python 3.7 or higher
* 4GB RAM minimum (8GB recommended)
* For large files, more memory may be required

Troubleshooting
--------------

Common Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Rust Compilation Fails**
   If the Rust implementation fails to compile, you can still use the pure Python implementation by setting ``use_rust=False`` in the API calls.

2. **Memory Issues**
   If you encounter memory issues during installation, try increasing your system's swap space or reducing the number of parallel compilation jobs.

3. **Dependency Conflicts**
   If you encounter dependency conflicts, consider using a virtual environment:

   .. code-block:: bash

       python -m venv hammock-env
       source hammock-env/bin/activate  # On Unix/macOS
       # or
       .\hammock-env\Scripts\activate  # On Windows
       pip install hammock 