#!/usr/bin/env python3
"""
Build the Rust HyperLogLog extension using maturin.
"""
import subprocess
import sys
import os
import shutil

def install_maturin():
    """Install maturin if not already installed."""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "maturin"])
        print("Maturin installed successfully")
    except subprocess.CalledProcessError:
        print("Failed to install maturin")
        sys.exit(1)

def build_extension():
    """Build the Rust extension using maturin."""
    os.chdir("rust_hll")
    try:
        # Build the extension in development mode
        subprocess.check_call([sys.executable, "-m", "maturin", "develop"])
        print("Rust extension built successfully")
    except subprocess.CalledProcessError:
        print("Failed to build Rust extension")
        sys.exit(1)

def verify_extension():
    """Verify the Rust extension was built correctly."""
    try:
        import rust_hll
        print("Successfully imported rust_hll module")
        print(f"Available objects: {dir(rust_hll)}")
        return True
    except ImportError as e:
        print(f"Failed to import rust_hll module: {e}")
        return False

if __name__ == "__main__":
    print("Installing maturin...")
    install_maturin()
    
    print("\nBuilding Rust extension...")
    build_extension()
    
    print("\nVerifying extension...")
    if verify_extension():
        print("\nBuild successful! You can now use FastHyperLogLog with Rust acceleration.")
    else:
        print("\nBuild failed. Using Python fallback.") 