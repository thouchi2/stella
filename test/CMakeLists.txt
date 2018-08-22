include_directories(../tools/single-block)

add_par_unit(test-models ${CMAKE_CURRENT_SOURCE_DIR}/test_models.cc ${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.cc)
target_link_libraries(test-models single-block stella)
set_par_sizes(test-models 1 3 4 8)

add_par_unit(test-models-3 ${CMAKE_CURRENT_SOURCE_DIR}/test_models_3.cc ${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.cc)
target_link_libraries(test-models-3 single-block stella)
set_par_sizes(test-models-3 1 3 4 8 9)

add_par_unit(test-signals ${CMAKE_CURRENT_SOURCE_DIR}/test_signals.cc ${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.cc)
target_link_libraries(test-signals single-block stella)
set_par_sizes(test-signals 1 3 4 8 10)

configure_file(../tools/single-block/PETScOptions.txt ${CMAKE_BINARY_DIR}/test/PETScOptions.txt COPYONLY)
