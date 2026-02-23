cmake_minimum_required(VERSION 3.16)
project(reed-solomon LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GMP
find_path(GMP_INCLUDE_DIR gmpxx.h
    HINTS /opt/homebrew/include /usr/local/include)
find_library(GMP_LIB gmp
    HINTS /opt/homebrew/lib /usr/local/lib)
find_library(GMPXX_LIB gmpxx
    HINTS /opt/homebrew/lib /usr/local/lib)

# Core library
add_library(rs_lib
    src/fp.cpp
    src/poly.cpp
)
target_include_directories(rs_lib PUBLIC src ${GMP_INCLUDE_DIR})
target_link_libraries(rs_lib PUBLIC ${GMP_LIB} ${GMPXX_LIB})

# Demo
add_executable(demo src/main.cpp)
target_link_libraries(demo PRIVATE rs_lib)

# Tests
add_executable(test src/test.cpp)
target_link_libraries(test PRIVATE rs_lib)
