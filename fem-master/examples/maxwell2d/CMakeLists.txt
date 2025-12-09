cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(example CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../example.cmake)

add_executable(example main.cpp ${EXTRA_INCS})
target_link_libraries(example ${EXTRA_LIBS})
