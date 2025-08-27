#!/usr/bin/env python3
"""
Build script for the C++ HyperLogLog extension.

This script helps build the C++ extension with proper error handling
and provides clear instructions for troubleshooting.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

def check_dependencies():
    """Check if required dependencies are available."""
    print("Checking dependencies...")
    
    # Check for Cython
    try:
        import Cython
        print(f"✓ Cython {Cython.__version__} found")
    except ImportError:
        print("✗ Cython not found")
        print("  Install with: pip install cython")
        return False
    
    # Check for numpy
    try:
        import numpy
        print(f"✓ NumPy {numpy.__version__} found")
    except ImportError:
        print("✗ NumPy not found")
        print("  Install with: pip install numpy")
        return False
    
    # Check for C++ compiler
    try:
        result = subprocess.run(['g++', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            print(f"✓ C++ compiler found: {version_line}")
        else:
            print("✗ C++ compiler (g++) not found or not working")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("✗ C++ compiler (g++) not found")
        print("  Install g++ or set up your C++ development environment")
        return False
    
    # Check for C++ source files
    hll_dir = Path("hll")
    required_files = ["hll.h", "hll.cpp", "kthread.h", "kthread.c"]
    
    for file in required_files:
        if not (hll_dir / file).exists():
            print(f"✗ Required C++ file not found: hll/{file}")
            return False
    
    print("✓ All C++ source files found")
    return True

def build_cpp_library():
    """Build the C++ library first."""
    print("\nBuilding C++ library...")
    
    hll_dir = Path("hll")
    os.chdir(hll_dir)
    
    try:
        # Clean previous builds
        subprocess.run(['make', 'clean'], check=True)
        
        # Build the library
        subprocess.run(['make', 'libhll.a'], check=True)
        
        # Test the build
        subprocess.run(['make', 'test'], check=True)
        
        print("✓ C++ library built successfully")
        os.chdir("..")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to build C++ library: {e}")
        os.chdir("..")
        return False

def build_python_extension():
    """Build the Python Cython extension."""
    print("\nBuilding Python Cython extension...")
    
    try:
        # Build the extension
        result = subprocess.run([
            sys.executable, 'setup.py', 'build_ext', '--inplace'
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Python extension built successfully")
            return True
        else:
            print(f"✗ Failed to build Python extension:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"✗ Error building Python extension: {e}")
        return False

def test_extension():
    """Test if the extension can be imported."""
    print("\nTesting extension import...")
    
    try:
        # Try to import the extension
        sys.path.insert(0, os.getcwd())
        from hammock.lib.cpp_hll_wrapper import CppHyperLogLog
        print("✓ C++ HyperLogLog extension imported successfully")
        
        # Test basic functionality
        hll = CppHyperLogLog(precision=12, hash_size=64)
        hll.add_string("test")
        cardinality = hll.cardinality()
        print(f"✓ Basic functionality test passed (cardinality: {cardinality:.2f})")
        
        return True
        
    except ImportError as e:
        print(f"✗ Failed to import extension: {e}")
        return False
    except Exception as e:
        print(f"✗ Extension test failed: {e}")
        return False

def install_extension():
    """Install the extension in development mode."""
    print("\nInstalling extension in development mode...")
    
    try:
        result = subprocess.run([
            sys.executable, '-m', 'pip', 'install', '-e', '.'
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Extension installed successfully")
            return True
        else:
            print(f"✗ Failed to install extension:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"✗ Error installing extension: {e}")
        return False

def main():
    """Main build process."""
    print("C++ HyperLogLog Extension Build Script")
    print("=" * 50)
    
    # Check dependencies
    if not check_dependencies():
        print("\nPlease install missing dependencies and try again.")
        return False
    
    # Build C++ library
    if not build_cpp_library():
        print("\nC++ library build failed. Please check the error messages above.")
        return False
    
    # Build Python extension
    if not build_python_extension():
        print("\nPython extension build failed. Please check the error messages above.")
        return False
    
    # Test extension
    if not test_extension():
        print("\nExtension test failed. Please check the error messages above.")
        return False
    
    # Install extension
    if not install_extension():
        print("\nExtension installation failed. Please check the error messages above.")
        return False
    
    print("\n" + "=" * 50)
    print("✓ C++ HyperLogLog extension built and installed successfully!")
    print("\nYou can now use the C++ HyperLogLog implementation in hammock:")
    print("  from hammock.lib.cpp_hyperloglog import CppHyperLogLogSketch")
    print("  sketch = CppHyperLogLogSketch(precision=12, hash_size=64)")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
