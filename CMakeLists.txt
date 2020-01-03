cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
project(hgauss LANGUAGES C VERSION 0.1.0)
set(PROJECT_DESCRIPTION "Implements the logarithm of the sum of exponentiations.")

# Generate compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Set RPATH only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

add_library(hgauss src/hgauss.c)
add_library(HGAUSS::hgauss ALIAS hgauss)

find_package(hcephes REQUIRED)

target_include_directories(hgauss
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(hgauss PRIVATE HCEPHES::hcephes)
set_target_properties(hgauss PROPERTIES PUBLIC_HEADER include/hgauss.h)
set_target_properties(hgauss PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(hgauss PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})
set_target_properties(hgauss PROPERTIES C_STANDARD 99)
set_target_properties(hgauss PROPERTIES DEFINE_SYMBOL HGAUSS_EXPORTS)
target_compile_definitions(hgauss PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>)
target_compile_definitions(hgauss PRIVATE $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)
target_link_libraries(hgauss INTERFACE $<$<BOOL:${UNIX}>:m>)
target_compile_options(hgauss PRIVATE
    $<$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>:/W3 /WX>
    $<$<NOT:$<STREQUAL:${CMAKE_C_COMPILER_ID},MSVC>>:-Wall -Wextra -Wstrict-prototypes
        -Wshadow -Wconversion -Wmissing-prototypes -Wno-unused-parameter>
)

install(TARGETS hgauss EXPORT hgauss-targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)
install(EXPORT hgauss-targets
    FILE hgauss-targets.cmake
    NAMESPACE HGAUSS::
    DESTINATION lib/cmake/hgauss
)
include(CMakePackageConfigHelpers)
configure_package_config_file(hgauss-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hgauss-config.cmake
    INSTALL_DESTINATION lib/cmake/hgauss
    NO_SET_AND_CHECK_MACRO
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hgauss-config.cmake
    DESTINATION lib/cmake/hgauss
)

enable_testing()
add_subdirectory(test)
