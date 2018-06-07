add_par_unit(test-models_3-3d ${CMAKE_CURRENT_SOURCE_DIR}/test_models_3.cc ${CMAKE_CURRENT_SOURCE_DIR}/../mpi_main.cc)
target_link_libraries(test-models_3-3d stella)
set_par_sizes(test-models_3-3d 4)
