cmake_minimum_required(VERSION 3.16)
project(chirp_engine VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_program(MKOCTFILE_EXECUTABLE mkoctfile)
if(NOT MKOCTFILE_EXECUTABLE)
    message(WARNING "mkoctfile not found in PATH. MEX compilation targets disabled.")
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/mex)

set(MEX_SOURCES
    src/mex/chirp_cfar_2d_mex.c
    src/mex/chirp_range_fft_mex.c
    src/mex/chirp_backprojection_mex.c
)

add_custom_target(chirp_mex ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Building MEX acceleration kernels via Makefile/mkoctfile"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

enable_testing()
add_test(NAME ChirpUnitTests
    COMMAND octave-cli --no-gui --eval "addpath('src'); addpath('tests'); run_all_tests;"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
