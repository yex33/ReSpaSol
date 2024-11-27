cmake_minimum_required(VERSION 3.1...3.28)

set(CMAKE_CXX_STANDARD 11)

add_library(mm_io mm_io.cpp mm_io.h)

add_library(libloadmatrix loadMatrixMarket.cpp loadMatrixMarket.h)
target_link_libraries(libloadmatrix PUBLIC mm_io)

