cmake_minimum_required(VERSION 3.20)
project(znmno2_model VERSION 0.1.0 LANGUAGES CXX Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The Fortran program uses bandsolver, fetched from GitHub by default;
# use a local checkout with -DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver
include(FetchContent)
FetchContent_Declare(bandsolver
  GIT_REPOSITORY https://github.com/jcbernard87/bandsolver.git
  GIT_TAG v0.1.2)
set(BANDSOLVER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BANDSOLVER_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BANDSOLVER_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(bandsolver)
# No fused multiply-adds anywhere, bandsolver included: the faithful modes reproduce the original programs exactly
if(TARGET bandsolver_f)
  target_compile_options(bandsolver_f PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffp-contract=off>)
endif()

add_subdirectory(fortran)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/CMakeLists.txt)
  add_subdirectory(cpp)
endif()
