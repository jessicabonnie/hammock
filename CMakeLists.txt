cmake_minimum_required(VERSION 3.18)
project(hammock LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(HAMMOCK_BUILD_TESTS "Build C++ unit tests" OFF)
option(HAMMOCK_BUILD_CLI "Build standalone hammock-cpp benchmark binary" ON)

# OpenMP for parallel mode-B/C processing.
find_package(OpenMP REQUIRED)

# pybind11 is supplied by scikit-build-core via the build-system requires
# in pyproject.toml. When invoked outside that pipeline (plain cmake), the
# user can still point at a system pybind11.
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found via find_package; the Python extension will be skipped. "
                   "Run via `pip install -e .` to build it.")
endif()

add_subdirectory(extern/hll)

set(HAMMOCK_CORE_SOURCES
    cpp/src/bagminhash_sketch.cpp
    cpp/src/bed_parser.cpp
    cpp/src/hll_sketch.cpp
    cpp/src/processing_modes.cpp
    cpp/src/stride.cpp
    cpp/src/xxhash.cpp
)

add_library(hammock_core STATIC ${HAMMOCK_CORE_SOURCES} $<TARGET_OBJECTS:hll_obj>)
target_include_directories(hammock_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)
target_link_libraries(hammock_core PUBLIC OpenMP::OpenMP_CXX hll_obj)
target_compile_options(hammock_core PRIVATE -O3 -march=native)

if(pybind11_FOUND)
    pybind11_add_module(_core MODULE bindings/_core.cpp)
    target_link_libraries(_core PRIVATE hammock_core)
    target_compile_options(_core PRIVATE -O3 -march=native)
    install(TARGETS _core LIBRARY DESTINATION hammock)
endif()

if(HAMMOCK_BUILD_CLI)
    add_executable(hammock-cpp cpp/app/hammock_cli.cpp)
    target_link_libraries(hammock-cpp PRIVATE hammock_core OpenMP::OpenMP_CXX)
    target_compile_options(hammock-cpp PRIVATE -O3 -march=native)
    # Not installed into the wheel; users who want it on $PATH can
    # `cmake --install` separately.
    install(TARGETS hammock-cpp RUNTIME DESTINATION bin OPTIONAL)
endif()

if(HAMMOCK_BUILD_TESTS)
    enable_testing()
    add_subdirectory(cpp/tests)
endif()
