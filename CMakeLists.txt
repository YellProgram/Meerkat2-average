cmake_minimum_required(VERSION 3.10)
project(Meerkat2-average)

#set(LIBHDF5_DIR lib/hdf5)
#
#set(HDF5_LIBRARIES hdf5_hl_cpp hdf5_cpp hdf5_hl hdf5 szip z )

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
    # Define a variable that controls the build architecture(s)
    # The user should set this via -DTARGET_ARCH=arm64 or -DTARGET_ARCH=x86_64
    set(TARGET_ARCH "arm64" CACHE STRING "Target architecture for macOS (e.g., arm64, x86_64).")

    # Set the main project's CMAKE_OSX_ARCHITECTURES based on the user's choice
    set(CMAKE_OSX_ARCHITECTURES "${TARGET_ARCH}" CACHE STRING "OSX Architectures" FORCE)

    # Define the argument string to pass to ExternalProject_Add
    set(ARCH_CMAKE_ARG "-DCMAKE_OSX_ARCHITECTURES=${TARGET_ARCH}")
else()
    set(ARCH_CMAKE_ARG "")
endif()

include(cmake/hdf5.cmake)

set (CMAKE_CXX_STANDARD 17)

set(KISSFFT_DIR lib/kiss_fft)
set(EIGEN_DIR lib/eigen-3.4.0)
include_directories(${EIGEN_DIR} ${KISSFFT_DIR} ${KISSFFT_DIR}/tools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") #-mtune=native -O1

set(SOURCE_FILES
        InputParameters.cpp InputParameters.h Exceptions.h Symmetry.cpp Symmetry.h IntensityData.cpp IntensityData.h)
add_library(meerkat-average-lib ${SOURCE_FILES})
DownloadHDF5For(meerkat-average-lib)

add_library(kissfft ${KISSFFT_DIR}/kiss_fft.c ${KISSFFT_DIR}/tools/kiss_fftnd.c) #
add_executable(Meerkat2-average main.cpp)
target_link_libraries(Meerkat2-average meerkat-average-lib  kissfft)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/lib/cxxtest/build_tools/cmake")
find_package(CxxTest)

if(CXXTEST_FOUND)
    include_directories(${CXXTEST_INCLUDES})
    cxx_test(test-meerkat-average test_meerkat_average.h)
    target_link_libraries(test-meerkat-average meerkat-average-lib)
endif()
