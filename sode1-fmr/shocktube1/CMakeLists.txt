cmake_minimum_required(VERSION 2.8)
if(POLICY CMP0037)
  cmake_policy(SET CMP0037 OLD)
endif()

project (formura_sin3)

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O2 -Wall -Wno-unknown-pragmas -mcmodel=large -std=gnu11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wno-unknown-pragmas -mcmodel=large")

add_custom_command(
  OUTPUT sode1-1.c
  COMMAND formura sode1-1.fmr --sleep 1
  COMMAND sleep 1
  DEPENDS sode1-1.fmr
)

add_library(shocktube shocktube.c)
include_directories("${PROJECT_SOURCE_DIR}")
FILE(GLOB LibSources sode1-1_internal*.c)
add_executable(main ${LibSources} sode1-1.c sode1-1-main.cpp)
target_link_libraries(main shocktube)
