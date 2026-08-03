cmake_minimum_required(VERSION 3.16)
project(highs_fortran_bridge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(HIGHS_ROOT "" CACHE PATH "Installed HiGHS prefix")
find_path(HIGHS_INCLUDE_DIR Highs.h
  PATHS "${HIGHS_ROOT}/include/highs" "${HIGHS_ROOT}/include")
find_library(HIGHS_LIBRARY NAMES highs libhighs
  PATHS "${HIGHS_ROOT}/lib" "${HIGHS_ROOT}/lib64")
if(NOT HIGHS_INCLUDE_DIR OR NOT HIGHS_LIBRARY)
  message(FATAL_ERROR "HiGHS headers/library not found under HIGHS_ROOT=${HIGHS_ROOT}")
endif()

add_library(highs_fortran_bridge SHARED highs_fortran_bridge.cpp)
target_include_directories(highs_fortran_bridge PRIVATE "${HIGHS_INCLUDE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/compat")
target_link_libraries(highs_fortran_bridge PRIVATE "${HIGHS_LIBRARY}")
set_target_properties(highs_fortran_bridge PROPERTIES OUTPUT_NAME "highs_fortran_bridge")

if(MINGW)
  target_link_options(highs_fortran_bridge PRIVATE -static-libgcc -static-libstdc++)
endif()
