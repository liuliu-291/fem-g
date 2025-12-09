## GmshFEM - Copyright (C) 2019-2020 Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# if CMAKE_BUILD_TYPE is specified use it; otherwise set the default build type
# to "Release" prior to calling project()
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type")
endif()

if(ENABLE_BESSEL)
  project(gmshfem CXX Fortran)
else()
  project(gmshfem CXX)
endif()

include(CTest)
enable_testing()

set(GMSHFEM_MAJOR_VERSION 1)
set(GMSHFEM_MINOR_VERSION 0)
set(GMSHFEM_PATCH_VERSION 0)
set(GMSHFEM_VERSION ${GMSHFEM_MAJOR_VERSION}.${GMSHFEM_MINOR_VERSION}.${GMSHFEM_PATCH_VERSION})

#----------------------------------
#
# S Y S T E M   C H E C K S
#
#----------------------------------

string(TIMESTAMP DATE "%Y%m%d")
if(NOT DATE)
  set(DATE "unknown")
endif()
set(GMSHFEM_DATE "${DATE}")

execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT HOSTNAME)
  set(HOSTNAME "unknown")
endif()
set(GMSHFEM_HOST "${HOSTNAME}")

execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT PACKAGER)
  set(PACKAGER "unknown")
endif()
string(REPLACE "\\" " " PACKAGER ${PACKAGER})
set(GMSHFEM_PACKAGER "${PACKAGER}")

if(APPLE)
  set(GMSHFEM_OS "macOS")
elseif(CYGWIN OR MSYS)
  # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
  # as a pure Windows build
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
     CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                    OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CXX_COMPILER_MACHINE MATCHES "mingw")
      set(GMSHFEM_OS "Windows")
      set(HAVE_WIN32 TRUE)
      add_definitions(-DWIN32)
      set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".lib" ".LIB" ".dll" ".DLL")
    endif()
  endif()
elseif(WIN32)
  set(GMSHFEM_OS "Windows")
  set(HAVE_WIN32 TRUE)
else()
  set(GMSHFEM_OS "${CMAKE_SYSTEM_NAME}")
endif()

#----------------------------------
#
# M A C R O S
#
#----------------------------------

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro()

macro(append_gmshfem_src DIRNAME FILES)
  set(LIST "")
  foreach(FILE ${FILES})
    list(APPEND LIST src/${DIRNAME}/${FILE})
  endforeach()
  set(GMSHFEM_SRC ${GMSHFEM_SRC};${LIST} PARENT_SCOPE)
  set(GMSHFEM_DIRS ${GMSHFEM_DIRS};${DIRNAME} PARENT_SCOPE)
endmacro()

macro(append_gmshfem_contrib DIRNAME FILES)
  set(LIST "")
  foreach(FILE ${FILES})
    list(APPEND LIST contrib/${DIRNAME}/${FILE})
  endforeach()
  set(GMSHFEM_SRC_CONTRIB ${GMSHFEM_SRC_CONTRIB};${LIST} PARENT_SCOPE)
  set(GMSHFEM_DIRS_CONTRIB ${GMSHFEM_DIRS_CONTRIB};${DIRNAME} PARENT_SCOPE)
endmacro()

macro(append_gmshfem_test FILES)
  set(LIST "")
  foreach(FILE ${FILES})
    list(APPEND LIST unitTests/coverage/${FILE})
  endforeach()
  set(GMSHFEM_SRC_TEST ${GMSHFEM_SRC_TEST};${LIST} PARENT_SCOPE)
endmacro()

macro(set_config_option STRING VARNAME)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "Found " ${STRING})
endmacro()

#----------------------------------
#
# O P T I O N S
#
#----------------------------------

set(CMAKE_CXX_STANDARD 17)

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-Wall" HAVE_WALL)
if(HAVE_WALL)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

# until Eigen is updated
check_cxx_compiler_flag("-Wno-unused-but-set-variable" HAVE_UNUSED)
if(HAVE_UNUSED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
endif()

# until Boost is updated and we remove all use of sprintf
check_cxx_compiler_flag("-Wno-deprecated-declarations" HAVE_DEPRECATED)
if(HAVE_DEPRECATED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

# remove "could not create compact unwind" linker warnings on macOS; but breaks
# exception handling
if(APPLE)
  # set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-no_compact_unwind")
endif()

opt(BESSEL "Enable Fortran Bessel functions" OFF)
opt(BOOST "Enable Boost library" ON)
opt(COVERAGE "Enable code coverage" OFF)
opt(EIGEN_USE_BLAS "Use BLAS with Eigen" OFF)
opt(MPI "Enable MPI" OFF)
opt(OPENMP "Enable OpenMP" ON)
opt(PETSC "Enable PETSc" ON)
opt(ROBINHOOD "Enable robin-hood hash-map" ON)
opt(SLEPC "Enable SLEPc" ON)
opt(SYSTEM_CONTRIB "Use system versions of contrib libraries, when possible" OFF)

#----------------------------------
#
# C O N S T A N T S
#
#----------------------------------

if(NOT GMSHFEM_DOF_FIELD_OFFSET)
  set(GMSHFEM_DOF_FIELD_OFFSET 1024)
endif()

if(NOT GMSHFEM_DOF_MEMORY_POOL_SIZE)
  set(GMSHFEM_DOF_MEMORY_POOL_SIZE 16384)
endif()

if(NOT GMSHFEM_FUNCTION_MEMORY_ALIGNMENT)
  set(GMSHFEM_FUNCTION_MEMORY_ALIGNMENT 512)
endif()

if(NOT GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE)
  set(GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE 1.5)
endif()

if(NOT GMSHFEM_FUNCTION_MEMORY_RESIDUAL)
  set(GMSHFEM_FUNCTION_MEMORY_RESIDUAL 1073741824) # 1GB
endif()

if(NOT GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE)
  set(GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE 16)
endif()

if(NOT GMSHFEM_NUMA_ALLOCATOR)
  set(GMSHFEM_NUMA_ALLOCATOR 0)
endif()

#----------------------------------
#
# A D D   D I R E C T O R I E S
#
#----------------------------------

add_subdirectory(src/algebra)
add_subdirectory(src/analytics)
add_subdirectory(src/common)
add_subdirectory(src/common/io)
add_subdirectory(src/equation)
add_subdirectory(src/dofs)
add_subdirectory(src/domain)
add_subdirectory(src/field)
add_subdirectory(src/function)
add_subdirectory(src/function/executionTree)
add_subdirectory(src/post)
add_subdirectory(src/post/mesh)
add_subdirectory(src/problem)
add_subdirectory(src/reorder)
add_subdirectory(src/system)
add_subdirectory(src/term)
add_subdirectory(src/term/evaluator)
include_directories(src/algebra src/analytics src/common src/common/io src/equation src/dofs src/domain src/field src/function src/function/executionTree src/post src/post/mesh src/problem src/reorder src/system src/term src/term/evaluator)

add_subdirectory(examples)
add_subdirectory(tutorials)

#----------------------------------
#
# O T H E R   L I B S
#
#----------------------------------

# Eigen (mandatory)

if(ENABLE_SYSTEM_CONTRIB)
  find_path(EIGEN_INC "Eigen/Dense" PATH_SUFFIXES eigen3)
  if(EIGEN_INC)
    include_directories(${EIGEN_INC})
    message(STATUS "Found Eigen")
    set(HAVE_CONTRIB_EIGEN FALSE)
  endif()
endif()
if(NOT EIGEN_INC)
  include_directories(contrib/eigen)
  message(STATUS "Found Eigen[contrib]")
  set(HAVE_CONTRIB_EIGEN TRUE)
endif()

# Gmsh (mandatory)

find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find Gmsh library")
endif()

find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find Gmsh header 'gmsh.h'")
endif()

set_config_option("Gmsh" HAVE_GMSH)

list(APPEND EXTRA_INCS ${GMSH_INC})
list(APPEND EXTRA_LIBS ${GMSH_LIB})

# PETSc

if(ENABLE_PETSC)

  if(PETSC_DIR)
    set(ENV_PETSC_DIR ${PETSC_DIR})
  else()
    set(ENV_PETSC_DIR $ENV{PETSC_DIR})
  endif()

  if(PETSC_ARCH)
    set(ENV_PETSC_ARCH ${PETSC_ARCH})
  else()
    set(ENV_PETSC_ARCH $ENV{PETSC_ARCH})
  endif()

  set(PETSC_POSSIBLE_CONF_FILES
      ${ENV_PETSC_DIR}/${ENV_PETSC_ARCH}/conf/petscvariables
      ${ENV_PETSC_DIR}/${ENV_PETSC_ARCH}/lib/petsc-conf/petscvariables
      ${ENV_PETSC_DIR}/${ENV_PETSC_ARCH}/lib/petsc/conf/petscvariables)

  foreach(FILE ${PETSC_POSSIBLE_CONF_FILES})
    if(EXISTS ${FILE})
      # old-style PETSc installations (using PETSC_DIR and PETSC_ARCH)
      message(STATUS "Using PETSc dir: ${ENV_PETSC_DIR}")
      message(STATUS "Using PETSc arch: ${ENV_PETSC_ARCH}")
      # find includes by parsing the petscvariables file
      file(STRINGS ${FILE} PETSC_VARIABLES NEWLINE_CONSUME)
    endif()
  endforeach()

  if(PETSC_VARIABLES)
    # try to find PETSC_CC_INCLUDES for PETSc >= 3.4
    string(REGEX MATCH "PETSC_CC_INCLUDES = [^\n\r]*" PETSC_PACKAGES_INCLUDES ${PETSC_VARIABLES})
    if(PETSC_PACKAGES_INCLUDES)
      string(REPLACE "PETSC_CC_INCLUDES = " "" PETSC_PACKAGES_INCLUDES ${PETSC_PACKAGES_INCLUDES})
    else()
      # try to find PETSC_PACKAGES_INCLUDES in older versions
      list(APPEND EXTRA_INCS ${ENV_PETSC_DIR}/include)
      list(APPEND EXTRA_INCS ${ENV_PETSC_DIR}/${ENV_PETSC_ARCH}/include)
      string(REGEX MATCH "PACKAGES_INCLUDES = [^\n\r]*" PETSC_PACKAGES_INCLUDES ${PETSC_VARIABLES})
      string(REPLACE "PACKAGES_INCLUDES = " "" PETSC_PACKAGES_INCLUDES ${PETSC_PACKAGES_INCLUDES})
    endif()

    if(PETSC_PACKAGES_INCLUDES)
      if(PETSC_PACKAGES_INCLUDES)
        string(REPLACE "-I" "" PETSC_PACKAGES_INCLUDES ${PETSC_PACKAGES_INCLUDES})
        string(REPLACE " " ";" PETSC_PACKAGES_INCLUDES ${PETSC_PACKAGES_INCLUDES})
        foreach(VAR ${PETSC_PACKAGES_INCLUDES})
          list(APPEND EXTRA_INCS ${VAR})
        endforeach()
      endif()
    endif()

    find_library(PETSC_LIBS petsc PATHS ${ENV_PETSC_DIR}/${ENV_PETSC_ARCH}/lib NO_DEFAULT_PATH)

    if(ENABLE_SLEPC)
      if(SLEPC_DIR)
        set(ENV_SLEPC_DIR ${SLEPC_DIR})
      else()
        set(ENV_SLEPC_DIR $ENV{SLEPC_DIR})
      endif()
      find_library(SLEPC_LIB slepc PATHS ${ENV_SLEPC_DIR}/${ENV_PETSC_ARCH}/lib NO_DEFAULT_PATH)
      if(SLEPC_LIB)
        find_path(SLEPC_INC "slepc.h" PATHS ${ENV_SLEPC_DIR} PATH_SUFFIXES include ${ENV_PETSC_ARCH}/include include/slepc NO_DEFAULT_PATH)
        if(SLEPC_INC)
          message(STATUS "Using SLEPc dir: ${ENV_SLEPC_DIR}")
          set_config_option("SLEPc" HAVE_SLEPC)
          list(APPEND EXTRA_LIBS ${SLEPC_LIB})
          list(APPEND EXTRA_INCS ${SLEPC_INC})
          find_path(SLEPC_INC2 "slepcconf.h" PATHS ${ENV_SLEPC_DIR} PATH_SUFFIXES ${ENV_PETSC_ARCH}/include NO_DEFAULT_PATH)
          if(SLEPC_INC2)
            list(APPEND EXTRA_INCS ${SLEPC_INC2})
          endif()
        endif()
      endif()
    endif()

    list(APPEND EXTRA_LIBS ${PETSC_LIBS})

    string(REGEX MATCH "PACKAGES_LIBS = [^\n\r]*" PLIBS ${PETSC_VARIABLES})
    if(PLIBS)
      string(REPLACE "PACKAGES_LIBS = " "" PLIBS ${PLIBS})
      string(STRIP ${PLIBS} PLIBS)
      list(APPEND EXTRA_LIBS "${PLIBS}")
    endif()

    string(REGEX MATCH "PETSC_EXTERNAL_LIB_BASIC = [^\n\r]*" PLIBS_BASIC ${PETSC_VARIABLES})
    if(PLIBS_BASIC)
      string(REPLACE "PETSC_EXTERNAL_LIB_BASIC = " "" PLIBS_BASIC ${PLIBS_BASIC})
      string(STRIP ${PLIBS_BASIC} PLIBS_BASIC)
      list(APPEND EXTRA_LIBS "${PLIBS_BASIC}")
    endif()

    string(REGEX MATCH "PCC_LINKER_LIBS = [^\n\r]*" LLIBS ${PETSC_VARIABLES})
    if(LLIBS)
      string(REPLACE "PCC_LINKER_LIBS = " "" LLIBS ${LLIBS})
      string(STRIP ${LLIBS} LLIBS)
      list(APPEND EXTRA_LIBS "${LLIBS}")
    endif()

    set_config_option("PETSc" HAVE_PETSC)
    list(APPEND EXTRA_LIBS ${PETSC_LIBS})
  else()
    find_library(PETSC_LIBS petsc)
    find_path(PETSC_INC "petsc.h" PATH_SUFFIXES include/petsc)
    if(PETSC_LIBS AND PETSC_INC)
      set_config_option("PETSc" HAVE_PETSC)
      list(APPEND EXTRA_LIBS ${PETSC_LIBS})
      list(APPEND EXTRA_INCS ${PETSC_INC})
    endif()
    if(ENABLE_SLEPC)
      find_library(SLEPC_LIBS slepc)
      find_path(SLEPC_INC "slepc.h" PATH_SUFFIXES include/slepc)
      if(SLEPC_LIBS AND SLEPC_INC)
        set_config_option("SLEPc" HAVE_SLEPC)
        list(APPEND EXTRA_LIBS ${SLEPC_LIBS})
        list(APPEND EXTRA_INCS ${SLEPC_INC})
      endif()
    endif()
  endif()

endif()

# Boost

if(ENABLE_BOOST)
  find_package(Boost)
  if(Boost_FOUND)
    list(APPEND EXTRA_INCS ${Boost_INCLUDE_DIRS})
    list(APPEND EXTRA_LIBS ${Boost_LIBRARIES})
    set_config_option("Boost" HAVE_BOOST)
  endif()
endif()

# Complex Bessel

if(ENABLE_BESSEL)
  add_subdirectory(contrib/bessel)
  include_directories(contrib/bessel)
  set_config_option("Bessel" HAVE_BESSEL)
endif()

# Eigen with BLAS

if(ENABLE_EIGEN_USE_BLAS)
  add_definitions(-DEIGEN_USE_BLAS)
  set_config_option("EigenWithBLAS" HAVE_EIGEN_USE_BLAS)
endif()

# MPI

if(ENABLE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND EXTRA_INCS ${MPI_CXX_INCLUDE_PATH})
    list(APPEND EXTRA_LIBS ${MPI_CXX_LIBRARIES})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    if(ENABLE_BESSEL)
      set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    endif()
    set_config_option("MPI" HAVE_MPI)
  endif()
endif()

# OpenMP

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()

  if(APPLE AND EXISTS "/usr/local/opt/libomp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    list(APPEND EXTRA_LIBS "-L/usr/local/opt/libomp/lib -lomp")
    set_config_option("OpenMP[Homebrew]" HAVE_OPENMP)
  elseif(APPLE AND EXISTS "/opt/local/lib/libomp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/local/include/libomp")
    list(APPEND EXTRA_LIBS "-L/opt/local/lib/libomp -lomp")
    set_config_option("OpenMP[MacPorts]" HAVE_OPENMP)
  elseif(APPLE AND EXISTS "/opt/homebrew/opt/libomp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    list(APPEND EXTRA_LIBS "-L/opt/homebrew/opt/libomp/lib -lomp")
    set_config_option("OpenMP[Homebrew]" HAVE_OPENMP)
  else()
    if(OPENMP_FOUND)
      set_config_option("OpenMP" HAVE_OPENMP)
    endif()
  endif()
endif()

# Robin-hood

if(ENABLE_ROBINHOOD)
  if(ENABLE_SYSTEM_CONTRIB)
    find_path(ROBINHOOD_INC "robin_hood.h" PATH_SUFFIXES robin-hood-hashing)
    if(ROBINHOOD_INC)
      include_directories(${ROBINHOOD_INC})
      set_config_option(HAVE_ROBINHOOD "RobinHood")
    endif()
  endif()
  if(NOT ROBINHOOD_INC)
    add_subdirectory(contrib/robin-hood-hashing)
    include_directories(contrib/robin-hood-hashing)
    set_config_option("RobinHood[contrib]" HAVE_ROBINHOOD)
  endif()
endif()

# Coverage

if(ENABLE_COVERAGE)
  set(CMAKE_BUILD_TYPE "Debug")
  include(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/coverage/CodeCoverage.cmake)
  append_coverage_compiler_flags()
  set(COVERAGE_EXCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/contrib/*"  "${CMAKE_CURRENT_SOURCE_DIR}/unitTests/coverage/*" "${CMAKE_CURRENT_SOURCE_DIR}/src/analytics/*" "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
  add_subdirectory(unitTests/coverage)
  include_directories(unitTests/coverage)
endif()

#----------------------------------
#
# F I N A L I Z E
#
#----------------------------------

list(SORT CONFIG_OPTIONS)
set(GMSHFEM_CONFIG_OPTIONS "")
foreach(OPT ${CONFIG_OPTIONS})
  set(GMSHFEM_CONFIG_OPTIONS "${GMSHFEM_CONFIG_OPTIONS} ${OPT}")
endforeach()

message(STATUS "")
message(STATUS "GmshFEM ${GMSHFEM_MAJOR_VERSION}.${GMSHFEM_MINOR_VERSION}.${GMSHFEM_PATCH_VERSION} has been configured with")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * Build options:" ${GMSHFEM_CONFIG_OPTIONS})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * C++ flags: " ${CMAKE_CXX_FLAGS})
if(ENABLE_BESSEL)
  message(STATUS " * Fortran compiler: " ${CMAKE_Fortran_COMPILER})
endif()
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS " * Constants: ")
message(STATUS "\t - GMSHFEM_DOF_FIELD_OFFSET = " ${GMSHFEM_DOF_FIELD_OFFSET})
message(STATUS "\t - GMSHFEM_DOF_MEMORY_POOL_SIZE = " ${GMSHFEM_DOF_MEMORY_POOL_SIZE})
message(STATUS "\t - GMSHFEM_FUNCTION_MEMORY_ALIGNMENT = " ${GMSHFEM_FUNCTION_MEMORY_ALIGNMENT})
message(STATUS "\t - GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE = " ${GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE})
message(STATUS "\t - GMSHFEM_FUNCTION_MEMORY_RESIDUAL = " ${GMSHFEM_FUNCTION_MEMORY_RESIDUAL})
message(STATUS "\t - GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE = " ${GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE})
message(STATUS "\t - GMSHFEM_NUMA_ALLOCATOR = " ${GMSHFEM_NUMA_ALLOCATOR})
message(STATUS "")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/common/gmshfemDefines.h.in ${CMAKE_CURRENT_BINARY_DIR}/common/gmshfemDefines.h)
include_directories(${EXTRA_INCS} ${CMAKE_CURRENT_BINARY_DIR}/common/)

# shared lib

include(GNUInstallDirs)
set(GMSHFEM_LIB ${CMAKE_INSTALL_LIBDIR})
set(GMSHFEM_INC ${CMAKE_INSTALL_INCLUDEDIR}/gmshfem)

if(NOT APPLE)
  set(CMAKE_SKIP_RPATH TRUE)
endif()

foreach(file ${GMSHFEM_SRC})
  string(REGEX MATCH ".+\\.h?h$" file ${file})
  if(file)
    list(APPEND HEADERS ${file})
  endif()
endforeach()
list(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/common/gmshfemDefines.h)

foreach(file ${GMSHFEM_SRC_CONTRIB})
  string(REGEX MATCH ".+\\.h?h$" file ${file})
  if(file)
    list(APPEND HEADERS_CONTRIB ${file})
  endif()
endforeach()

add_library(shared SHARED "${GMSHFEM_SRC};${GMSHFEM_SRC_CONTRIB}")
set_target_properties(shared PROPERTIES OUTPUT_NAME gmshfem VERSION ${GMSHFEM_MAJOR_VERSION}.${GMSHFEM_MINOR_VERSION}.${GMSHFEM_PATCH_VERSION} SOVERSION ${GMSHFEM_MAJOR_VERSION}.${GMSHFEM_MINOR_VERSION})

if(WIN32)
  set_target_properties(shared PROPERTIES PREFIX "" IMPORT_PREFIX "" IMPORT_SUFFIX ".lib" COMPILE_FLAGS "-DGMSHFEM_DLL -DGMSHFEM_DLL_EXPORT")
endif()
target_link_libraries(shared ${EXTRA_LIBS})

# Windows specific linker options
if(WIN32 AND NOT MSVC)
  set_target_properties(shared PROPERTIES LINK_FLAGS "-static")
endif()

install(TARGETS shared DESTINATION ${GMSHFEM_LIB})
install(FILES ${HEADERS} DESTINATION ${GMSHFEM_INC})
install(FILES ${HEADERS_CONTRIB} DESTINATION ${GMSHFEM_INC})
if(HAVE_CONTRIB_EIGEN)
  install(DIRECTORY contrib/eigen/Eigen DESTINATION ${GMSHFEM_INC})
endif()

# binary

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
  add_executable(gmshfem ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
  target_link_libraries(gmshfem shared)
endif()

# coverage

if(ENABLE_COVERAGE)
  add_executable(coverage ${GMSHFEM_SRC_TEST})
  target_link_libraries(coverage shared)

  setup_target_for_coverage_gcovr_html(NAME coverage_report EXECUTABLE ./coverage -debug DEPENDENCIES coverage)
endif()
