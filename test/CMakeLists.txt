include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_hgauss test_hgauss.c)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_hgauss hgauss)
else()
  target_link_libraries (test_hgauss hgauss m)
endif()
