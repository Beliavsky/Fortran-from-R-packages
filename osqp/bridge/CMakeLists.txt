cmake_minimum_required(VERSION 3.18)
project(osqp_fortran_bridge LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(OSQP_ROOT "" CACHE PATH "Installed OSQP prefix")
find_path(OSQP_INCLUDE_DIR osqp.h
  PATHS "${OSQP_ROOT}/include/osqp" "${OSQP_ROOT}/include")
find_library(OSQP_LIBRARY NAMES osqpstatic osqp libosqpstatic libosqp
  PATHS "${OSQP_ROOT}/lib" "${OSQP_ROOT}/lib64")
if(NOT OSQP_INCLUDE_DIR OR NOT OSQP_LIBRARY)
  message(FATAL_ERROR "OSQP headers/library not found under OSQP_ROOT=${OSQP_ROOT}")
endif()

add_library(osqp_fortran_bridge SHARED osqp_fortran_bridge.c)
target_include_directories(osqp_fortran_bridge PRIVATE
  "${OSQP_INCLUDE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(osqp_fortran_bridge PRIVATE "${OSQP_LIBRARY}")
set_target_properties(osqp_fortran_bridge PROPERTIES OUTPUT_NAME "osqp_fortran_bridge")

if(MINGW)
  target_link_options(osqp_fortran_bridge PRIVATE -static-libgcc)
endif()
