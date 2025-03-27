#!/usr/bin/env python
import os
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Define the build command for the Rust extension
class RustExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])

class RustBuildCommand(build_ext):
    def build_extension(self, ext):
        if not isinstance(ext, RustExtension):
            super().build_extension(ext)
            return
        
        # Build the Rust extension
        try:
            # Check if Rust and Cargo are installed
            subprocess.check_call(["cargo", "--version"])
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Rust or Cargo not found. Please install Rust: https://www.rust-lang.org/tools/install")
            sys.exit(1)
        
        # Set the output directory
        out_dir = os.path.dirname(os.path.abspath(self.get_ext_fullpath(ext.name)))
        os.makedirs(out_dir, exist_ok=True)
        
        # Build the extension
        target_dir = os.path.abspath(".")
        build_cmd = ["cargo", "build", "--release"]
        subprocess.check_call(build_cmd, cwd=target_dir)
        
        # Copy the compiled extension to the output directory
        ext_path = self.get_ext_fullpath(ext.name)
        ext_dir = os.path.dirname(ext_path)
        os.makedirs(ext_dir, exist_ok=True)
        
        # Find the shared library file
        if sys.platform == "darwin":
            lib_name = "librust_hll.dylib"
        elif sys.platform == "win32":
            lib_name = "rust_hll.dll"
        else:
            lib_name = "librust_hll.so"
        
        lib_path = os.path.join(target_dir, "target", "release", lib_name)
        
        # Rename the library to match what Python expects
        if sys.platform == "darwin":
            dest_name = os.path.join(ext_dir, "rust_hll.so")
        elif sys.platform == "win32":
            dest_name = os.path.join(ext_dir, "rust_hll.pyd")
        else:
            dest_name = os.path.join(ext_dir, "rust_hll.so")
        
        # Copy the file
        import shutil
        shutil.copyfile(lib_path, dest_name)
        print(f"Copied extension from {lib_path} to {dest_name}")

# Setup the extension
setup(
    name="rust_hll",
    version="0.1.0",
    ext_modules=[RustExtension("rust_hll")],
    cmdclass={"build_ext": RustBuildCommand},
)

if __name__ == "__main__":
    # If run directly, build the extension
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    subprocess.check_call([sys.executable, "build.py", "build_ext", "--inplace"])
    print("Rust HyperLogLog extension built successfully!") 