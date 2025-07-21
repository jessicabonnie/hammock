#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import os
import warnings

with open("README.md", "r") as fh:
    long_description = fh.read()

# Try to import Cython for optional acceleration
try:
    from Cython.Build import cythonize
    import numpy as np
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

class OptionalBuildExt(build_ext):
    """Build extension that handles Cython failures gracefully."""
    
    def run(self):
        try:
            super().run()
        except Exception as e:
            warnings.warn(
                f"Failed to build Cython extensions: {e}\n"
                f"Hammock will work but without Cython acceleration.\n"
                f"To enable acceleration, install Cython: pip install cython",
                UserWarning
            )

def get_extensions():
    """Get list of extensions to build."""
    if not USE_CYTHON:
        return []
    
    extensions = [
        Extension(
            "hammock.lib._hyperloglog_ext",
            sources=["hammock/lib/_hyperloglog_ext.pyx"],
            include_dirs=[np.get_include()],
            extra_compile_args=["-O3"],
            extra_link_args=[],
        )
    ]
    
    return cythonize(extensions, compiler_directives={
        'language_level': "3",
        'boundscheck': False,
        'wraparound': False,
        'cdivision': True,
        'embedsignature': True,
    })

# Automatically find all scripts in the scripts directory
def get_scripts():
    scripts_dir = 'scripts'
    if os.path.exists(scripts_dir):
        scripts = []
        for file in os.listdir(scripts_dir):
            if not file.startswith('.') and not file.startswith('__'):
                # Include all files except hidden files and Python cache
                file_path = os.path.join(scripts_dir, file)
                if os.path.isfile(file_path):
                    scripts.append(file_path)
        return scripts
    return []

setup(
    name="hammock",
    version="0.2.2",
    author="Jessica Bonnie",
    author_email="jbonnie@jhu.edu",
    description="Cardinality Estimation and Interval Sketches",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jessicabonnie/hammock",
    packages=find_packages(),
    ext_modules=get_extensions(),
    cmdclass={'build_ext': OptionalBuildExt},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=get_scripts(),  # Automatically include all scripts
    python_requires='>=3.6',
    install_requires=[
        "numpy", 
        "matplotlib",
        "scipy",
        "xxhash",  # Required for FastHyperLogLog
    ],
    extras_require={
        'fast': ['cython>=0.29.0'],  # Optional Cython acceleration
        'rust': ['maturin>=1.0.0'],
        'dev': [
            'pytest>=7.0.0',
            'pytest-benchmark>=4.0.0',
            'psutil>=5.0.0',  # Required for memory monitoring tests
            'black>=22.0.0',
            'isort>=5.0.0',
            'mypy>=0.900',
            'maturin>=1.0.0',
            'cython>=0.29.0',  # Include Cython in dev requirements
        ],
    },
    entry_points={
        'console_scripts': [
            'hammock=hammock.hammock:main',
        ],
    },
) 