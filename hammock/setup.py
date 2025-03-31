from setuptools import setup, find_packages

setup(
    name="hammock",
    version="0.1.0",
    packages=find_packages(),
    package_dir={"": "."},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=[
        "setuptools",
    ],
) 