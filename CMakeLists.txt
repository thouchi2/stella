cmake_minimum_required(VERSION 2.8)
project(stella)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(cxx11)
include(cxx_config)
include(testing)

option(ENABLE_UNIT_TESTS "Enable unit testing" OFF)
option(AutodetectMPI "Automatically detect MPI" ON)

enable_language(C CXX Fortran)

find_package(PETSc REQUIRED)

set(req-lib ${PETSC_LIBRARIES} m)
set(req-inc ${PETSC_INCLUDES})

if(AutodetectMPI)
  find_package(MPI REQUIRED)
  list(APPEND req-lib ${MPI_C_LIBRARIES})
  list(APPEND req-inc ${MPI_C_INCLUDE_PATH})
endif()

include_directories(src)

if(ENABLE_UNIT_TESTS)
  add_subdirectory(test)
endif()

add_subdirectory(src)
add_library(stella ${stella_sources})
target_compile_definitions(stella PUBLIC PETSC_OPTIONS_TAKES_DATABASE)
target_link_libraries(stella PUBLIC ${req-lib})
target_include_directories(stella PUBLIC ${req-inc})

add_subdirectory(tools)
