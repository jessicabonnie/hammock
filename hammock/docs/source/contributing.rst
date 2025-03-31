Contributing to Hammock
======================

Thank you for your interest in contributing to Hammock! This document provides guidelines and instructions for contributing to the project.

Development Setup
---------------

1. Fork the repository on GitHub
2. Clone your fork locally:

   .. code-block:: bash

       git clone https://github.com/YOUR_USERNAME/hammock.git
       cd hammock

3. Install development dependencies:

   .. code-block:: bash

       pip install -e ".[dev]"

4. Set up pre-commit hooks (optional but recommended):

   .. code-block:: bash

       pre-commit install

Code Style
----------

Hammock follows these coding standards:

* Python code should follow PEP 8 guidelines
* Use type hints for all function parameters and return values
* Document all public functions and classes using docstrings
* Use the Google style for docstrings

Example:

.. code-block:: python

    def compare_bed_files(file1: str, file2: str, **kwargs) -> Dict[str, float]:
        """Compare two BED files and return similarity metrics.

        Args:
            file1: Path to first BED file
            file2: Path to second BED file
            **kwargs: Additional arguments passed to compare_files()

        Returns:
            Dictionary containing similarity metrics
        """
        pass

Testing
-------

* Write tests for all new features
* Ensure all tests pass before submitting a pull request
* Run tests with:

  .. code-block:: bash

      pytest

Documentation
------------

* Update documentation for any new features
* Add examples for new functionality
* Ensure all public APIs are documented
* Build documentation locally:

  .. code-block:: bash

      cd docs
      make html

Submitting Changes
----------------

1. Create a new branch for your feature:

   .. code-block:: bash

       git checkout -b feature/your-feature-name

2. Make your changes and commit them:

   .. code-block:: bash

       git commit -m "Add your feature"

3. Push to your fork:

   .. code-block:: bash

       git push origin feature/your-feature-name

4. Create a pull request on GitHub

Pull Request Guidelines
---------------------

* Provide a clear description of your changes
* Include tests for new features
* Update documentation as needed
* Ensure all CI checks pass
* Request review from maintainers

Code Review Process
-----------------

1. Maintainers will review your pull request
2. Address any feedback and make requested changes
3. Once approved, your changes will be merged

Reporting Issues
--------------

When reporting issues, please include:

* A clear description of the problem
* Steps to reproduce the issue
* Expected behavior
* Actual behavior
* Your environment (Python version, OS, etc.)

Getting Help
-----------

* Open an issue on GitHub
* Join our community discussions
* Contact the maintainers

Thank you for contributing to Hammock! 