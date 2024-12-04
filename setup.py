from setuptools import setup, find_packages

setup(
    name="hammock",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "xxhash",
        "scipy",
        "memory_profiler",
        "pandas",
        "matplotlib",
    ],
    entry_points={
        'console_scripts': [
            'hammock=hammock.hammock:main',
        ],
    },
) 