cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(tuto CXX)

find_package(MPI)
if(MPI_FOUND) 
  list(APPEND EXTRA_INCS ${MPI_CXX_INCLUDE_PATH})
  list(APPEND EXTRA_LIBS ${MPI_CXX_LIBRARIES})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/../tutorial.cmake)

add_executable(tuto main.cpp ${EXTRA_INCS})
target_link_libraries(tuto ${EXTRA_LIBS})
