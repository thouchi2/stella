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
if(AutoDetectMPI)
  find_package(MPI REQUIRED)
endif()

include_directories(${MPI_INCLUDE_PATH} ${PETSC_INCLUDES} src)

if(ENABLE_UNIT_TESTS)
  add_subdirectory(test)
endif()

add_subdirectory(src)
add_library(stella ${stella_sources})
target_compile_definitions(stella PUBLIC PETSC_OPTIONS_TAKES_DATABASE)
target_link_libraries(stella PUBLIC
  ${MPI_C_LIBRARIES}
  ${PETSC_LIBRARIES}
  m
)

add_subdirectory(tools)
