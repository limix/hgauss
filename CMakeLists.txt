cmake_minimum_required (VERSION 2.6.4)

project(hgauss C)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)

  if (GCC_VERSION VERSION_GREATER 4.9)
    set (CMAKE_C_STANDARD 11)
  else ()
    set (CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")
  endif ()
endif ()

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-Dinline=__inline)
endif()

file(STRINGS "VERSION" HGAUSS_VERSION)

include_directories( include/ )
add_subdirectory(hgauss)

add_library(hgauss SHARED hgauss/hgauss.c)

add_library(hgauss_static STATIC hgauss/hgauss.c)

set_target_properties(hgauss
  PROPERTIES
    VERSION ${HGAUSS_VERSION}
PUBLIC_HEADER include/hgauss/hgauss.h)

install(TARGETS hgauss
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
    PUBLIC_HEADER DESTINATION include/hgauss)

install(TARGETS hgauss_static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/hgauss)

add_subdirectory(test)

enable_testing()
add_test(test_hgauss test/test_hgauss)
