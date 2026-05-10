# Vendored copy of mindis/hll (https://github.com/mindis/hll).
# Compiled as an OBJECT library so the resulting object files are linked
# directly into hammock_core.

add_library(hll_obj OBJECT
    hll.cpp
    kthread.c
)

target_include_directories(hll_obj PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# THREADSAFE selects the lock-protected add() path; required because the
# unguarded operator++ has races on the register-update step.
target_compile_definitions(hll_obj PUBLIC THREADSAFE)

set_source_files_properties(hll.cpp PROPERTIES COMPILE_FLAGS "-std=c++17")

target_compile_options(hll_obj PRIVATE -O3 -march=native)

# PIC so these objects can be linked into the pybind11 shared module.
set_target_properties(hll_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
