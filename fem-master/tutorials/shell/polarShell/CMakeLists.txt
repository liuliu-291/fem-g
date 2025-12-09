cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(tuto CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../tutorial.cmake)

add_executable(tuto main.cpp ${EXTRA_INCS})
target_link_libraries(tuto ${EXTRA_LIBS})
