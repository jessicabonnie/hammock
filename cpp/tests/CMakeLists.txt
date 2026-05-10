include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
)
FetchContent_MakeAvailable(doctest)

add_executable(hammock_tests
    test_main.cpp
    test_bed_parser.cpp
    test_hll_sketch.cpp
    test_stride.cpp
    test_modes.cpp
)
target_link_libraries(hammock_tests PRIVATE hammock_core doctest::doctest OpenMP::OpenMP_CXX)
target_compile_options(hammock_tests PRIVATE -O2)

add_test(NAME hammock_tests COMMAND hammock_tests)
