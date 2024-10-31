# Specify minimum required version of CMake
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Specify project, compiler and languages
project(assignment4 LANGUAGES CXX)

# Include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++20")

set(SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/discrete_cosine_transform.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stream.cpp
)

add_executable(uvid_compress ${CMAKE_CURRENT_SOURCE_DIR}/src/uvid_compress.cpp ${SOURCES})
add_executable(uvid_decompress ${CMAKE_CURRENT_SOURCE_DIR}/src/uvid_decompress.cpp ${SOURCES})
add_executable(huffman ${CMAKE_CURRENT_SOURCE_DIR}/src/huffman.cpp)