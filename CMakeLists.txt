cmake_minimum_required(VERSION 3.14)
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()
project(Meerkat2-average)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
    set(TARGET_ARCH "arm64" CACHE STRING "Target architecture for macOS (e.g., arm64, x86_64).")
    set(CMAKE_OSX_ARCHITECTURES "${TARGET_ARCH}" CACHE STRING "OSX Architectures" FORCE)

    # Argument that will be passed to hdf5 and bitshuffle
    set(ARCH_CMAKE_ARG "-DCMAKE_OSX_ARCHITECTURES=${TARGET_ARCH}")
else()
    set(ARCH_CMAKE_ARG "")
endif()

include(cmake/hdf5.cmake)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

# Eigen - header-only linear algebra library
FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
)
FetchContent_Populate(eigen)

# kiss_fft - small FFT library
FetchContent_Declare(
    kiss_fft
    GIT_REPOSITORY https://github.com/mborgerding/kissfft.git
    GIT_TAG v131
)
FetchContent_Populate(kiss_fft)

include_directories(${eigen_SOURCE_DIR} ${kiss_fft_SOURCE_DIR} ${kiss_fft_SOURCE_DIR}/tools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") #-mtune=native -O1

set(SOURCE_FILES
        InputParameters.cpp InputParameters.h Exceptions.h Symmetry.cpp Symmetry.h IntensityData.cpp IntensityData.h)
add_library(meerkat-average-lib ${SOURCE_FILES})
DownloadHDF5For(meerkat-average-lib)

add_library(kissfft ${kiss_fft_SOURCE_DIR}/kiss_fft.c ${kiss_fft_SOURCE_DIR}/tools/kiss_fftnd.c)
target_include_directories(kissfft PUBLIC ${kiss_fft_SOURCE_DIR} ${kiss_fft_SOURCE_DIR}/tools)

add_executable(Meerkat2-average main.cpp)
target_link_libraries(Meerkat2-average meerkat-average-lib kissfft)

# CxxTest
FetchContent_Declare(
    cxxtest
    GIT_REPOSITORY https://github.com/CxxTest/cxxtest.git
    GIT_TAG master
)
FetchContent_Populate(cxxtest)

set(CXXTEST_PATH ${cxxtest_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${cxxtest_SOURCE_DIR}/build_tools/cmake")
find_package(CxxTest)

if(CXXTEST_FOUND)
    include_directories(${CXXTEST_INCLUDES})
    cxx_test(test-meerkat-average test_meerkat_average.h)
    target_link_libraries(test-meerkat-average meerkat-average-lib)
endif()
