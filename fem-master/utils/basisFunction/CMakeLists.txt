cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(demo CXX)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../tutorials/tutorial.cmake)

add_executable(demo main.cpp ${EXTRA_INCS})
target_link_libraries(demo ${EXTRA_LIBS})
