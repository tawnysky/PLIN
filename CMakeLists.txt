cmake_minimum_required(VERSION 3.16)

project(Plin)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS  "-Iinclude -mclflushopt -mclwb -march=native -O3")

add_executable(${PROJECT_NAME} test.cpp)

target_link_libraries(${PROJECT_NAME} -lpmemobj -pthread)
