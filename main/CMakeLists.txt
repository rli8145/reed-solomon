cmake_minimum_required(VERSION 3.16)
project(basepoint LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GMP
find_path(GMP_INCLUDE_DIR gmpxx.h
    HINTS /opt/homebrew/include /usr/local/include)
find_library(GMP_LIB gmp
    HINTS /opt/homebrew/lib /usr/local/lib)
find_library(GMPXX_LIB gmpxx
    HINTS /opt/homebrew/lib /usr/local/lib)

# Library
add_library(cat_lib
    src/fp.cpp
)
target_include_directories(cat_lib PUBLIC src ${GMP_INCLUDE_DIR})
target_link_libraries(cat_lib PUBLIC ${GMP_LIB} ${GMPXX_LIB})

# Demo
add_executable(demo src/main.cpp)
target_link_libraries(demo PRIVATE cat_lib)

# Tests (Google Test)
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(test_fp src/test_fp.cpp)
target_link_libraries(test_fp PRIVATE cat_lib GTest::gtest_main)
add_test(NAME test_fp COMMAND test_fp)
