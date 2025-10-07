#!/bin/bash
# Quick script to build and run the example

set -e  # Exit on error

echo "=== Building Hammock Mode A ==="
make

echo ""
echo "=== Running Example ==="
cd examples
../bin/hammock_mode_a files_list.txt primary_list.txt -p 14 -o example_output.csv -v

echo ""
echo "=== Results ==="
cat example_output.csv

echo ""
echo "=== Done! ==="
echo "Binary location: bin/hammock_mode_a"
echo "Example output: examples/example_output.csv"

