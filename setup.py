#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hammock",
    version="0.2.1",
    author="Jessica Bonnie",
    author_email="jbonnie@jhu.edu",
    description="Cardinality Estimation and Interval Sketches",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jessicabonnie/hammock",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['scripts/encode_heatmap.R', 'scripts/compare_sim_matrices.py'],
    python_requires='>=3.6',
    install_requires=[
        "numpy", 
        "matplotlib",
        "scipy"
    ],
    extras_require={
        'rust': ['maturin>=1.0.0'],
        'dev': [
            'pytest>=7.0.0',
            'pytest-benchmark>=4.0.0',
            'black>=22.0.0',
            'isort>=5.0.0',
            'mypy>=0.900',
            'maturin>=1.0.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'hammock=hammock.hammock:main',
        ],
    },
) 