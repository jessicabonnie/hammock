from setuptools import setup, find_packages

setup(
    name="hammock",
    version="0.1.0",
    description="A library for sketching and comparing files containing lists of things, with special capabilities for comparing bedfiles",
    author="Jessica Bonnie",
    author_email="jbonnie1@jhu.edu",
    packages=find_packages(),
    package_dir={"": "."},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=[
        "setuptools",
        "numpy>=1.19.0",
        "pandas>=1.2.0",
        "scipy>=1.6.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "mypy>=0.800",
            "flake8>=3.9",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics genomics bedfiles sketching comparison",
    url="https://github.com/jbonnie1/hammock",
) 