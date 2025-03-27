from setuptools import setup, find_packages

setup(
    name="hammock",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "xxhash>=2.0.0",
        "pyBigWig>=0.3.18",
        "biopython>=1.79",
        "pysam>=0.22.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.0.0",
            "pytest-timeout>=2.1.0",
            "mypy>=0.900",
        ],
    },
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'hammock=hammock.hammock:main',
        ],
    },
) 