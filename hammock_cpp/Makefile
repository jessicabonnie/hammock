CXX=g++
CC=gcc

# Paths
HLL_DIR=../hll
SRC_DIR=src
BUILD_DIR=build
BIN_DIR=bin

# Compiler flags
FLAGS=-O3 -funroll-loops -pipe -march=native -I$(HLL_DIR) -std=c++17 -Wall -Wextra -Wdisabled-optimization -DNDEBUG -DTHREADSAFE -fopenmp

ifeq ($(shell uname),Darwin)
	FLAGS := $(FLAGS) -Wa,-q
endif

# Create directories
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

# Default target
all: hammock

# Build HLL library first if needed
$(HLL_DIR)/hll.o:
	$(MAKE) -C $(HLL_DIR) hll.o

# Build hammock (supports both Mode A and B)
hammock: $(HLL_DIR)/hll.o $(SRC_DIR)/hammock.cpp
	$(CXX) $(FLAGS) $(HLL_DIR)/hll.o $(SRC_DIR)/hammock.cpp -o $(BIN_DIR)/hammock

# Install target (copies binary to root of hammock_cpp)
install: hammock
	cp $(BIN_DIR)/hammock ./hammock

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/* hammock test_hammock

# Clean everything including HLL objects
clean-all: clean
	$(MAKE) -C $(HLL_DIR) clean

# Build and run tests (Modes A and B)
test: $(HLL_DIR)/hll.o test/test_hammock.cpp
	$(CXX) $(FLAGS) -DHAMMOCK_TEST_MODE $(HLL_DIR)/hll.o test/test_hammock.cpp -o $(BIN_DIR)/test_hammock
	$(BIN_DIR)/test_hammock

# Run example with Mode A
example_a: hammock
	@echo "=== Mode A Example (Interval-based) ==="
	cd examples && ../$(BIN_DIR)/hammock files_list.txt primary_list.txt --mode A -p 14 -o example_output_A.csv -v

# Run example with Mode B
example_b: hammock
	@echo "=== Mode B Example (Point-based) ==="
	cd examples && ../$(BIN_DIR)/hammock files_list.txt primary_list.txt --mode B -p 14 -o example_output_B.csv -v

# Run example with Mode C
example_c: hammock
	@echo "=== Mode C Example (Combined) ==="
	cd examples && ../$(BIN_DIR)/hammock files_list.txt primary_list.txt --mode C -p 14 -o example_output_C.csv -v

# Run all examples
example: example_a example_b example_c

.PHONY: all clean clean-all install example example_a example_b example_c test

