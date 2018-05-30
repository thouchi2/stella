set(single-block_sources
  ${CMAKE_CURRENT_SOURCE_DIR}/grid.c
  ${CMAKE_CURRENT_SOURCE_DIR}/map.c
  ${CMAKE_CURRENT_SOURCE_DIR}/option.c
  ${CMAKE_CURRENT_SOURCE_DIR}/problem.c
  ${CMAKE_CURRENT_SOURCE_DIR}/solver.c
  ${CMAKE_CURRENT_SOURCE_DIR}/state.c
  ${CMAKE_CURRENT_SOURCE_DIR}/boundary.c
)

include_directories(.)
add_library(single-block ${single-block_sources})
target_link_libraries(single-block stella)

configure_file(PETScOptions.txt ${CMAKE_BINARY_DIR}/PETScOptions.txt COPYONLY)
