add_subdirectory(single-block)

add_executable(driver-sb ${CMAKE_CURRENT_SOURCE_DIR}/single-block/driver.c)
target_link_libraries(driver-sb single-block stella)
