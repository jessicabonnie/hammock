#!/usr/bin/env python3
import os
import sys

def debug_rust_path():
    """Debug the path to the Rust extension."""
    # Get current file path
    current_file = os.path.abspath(__file__)
    print(f"Current file: {current_file}")
    
    # Get directory of current file
    current_dir = os.path.dirname(current_file)
    print(f"Current directory: {current_dir}")
    
    # Path calculation for rust_hll directory
    rust_hll_dir = os.path.join(current_dir, "rust_hll")
    print(f"Expected rust_hll directory: {rust_hll_dir}")
    print(f"Directory exists: {os.path.exists(rust_hll_dir)}")
    
    # Check for Rust library file
    lib_files = []
    if os.path.exists(rust_hll_dir):
        target_dir = os.path.join(rust_hll_dir, "target", "release")
        if os.path.exists(target_dir):
            for filename in os.listdir(target_dir):
                if "rust_hll" in filename and (filename.endswith(".so") or filename.endswith(".dylib") or filename.endswith(".dll")):
                    lib_files.append(os.path.join(target_dir, filename))
    
    print(f"Found library files: {lib_files}")
    
    # Try to add rust_hll to path and import
    sys.path.append(rust_hll_dir)
    try:
        import rust_hll
        print(f"Successfully imported rust_hll")
        print(f"Module path: {rust_hll.__file__}")
    except ImportError as e:
        print(f"Failed to import rust_hll: {e}")
    
    # Try alternative approach - look for compiled extension in target directory
    if lib_files:
        for lib_file in lib_files:
            target_dir = os.path.dirname(lib_file)
            if target_dir not in sys.path:
                sys.path.append(target_dir)
        
        try:
            import rust_hll
            print(f"Successfully imported rust_hll from target directory")
            print(f"Module path: {rust_hll.__file__}")
        except ImportError as e:
            print(f"Failed to import rust_hll from target directory: {e}")
    
    # Print Python path
    print("\nPython path:")
    for path in sys.path:
        print(f"  {path}")

if __name__ == "__main__":
    debug_rust_path() 