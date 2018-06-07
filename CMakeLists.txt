cmake_minimum_required(VERSION 2.8)
project(stella)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(cxx11)
include(cxx_config)
include(testing)

enable_language(C CXX Fortran)

find_package(PETSc REQUIRED)

include_directories(${MPI_INCLUDE_PATH} ${PETSC_INCLUDES} src)

add_subdirectory(src)
add_library(stella ${stella_sources})
target_compile_definitions(stella PUBLIC PETSC_OPTIONS_TAKES_DATABASE)
target_link_libraries(stella PUBLIC
  ${MPI_C_LIBRARIES}
  ${PETSC_LIBRARIES}
  m
)

add_subdirectory(tools)
