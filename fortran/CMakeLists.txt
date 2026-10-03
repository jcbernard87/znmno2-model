add_executable(znmno2_f params.f90 equilibria.f90 corrected.f90 faithful.f90 main.f90)
target_link_libraries(znmno2_f PRIVATE bandsolver_f)
target_include_directories(znmno2_f PRIVATE ${CMAKE_BINARY_DIR}/modules)
# No fused multiply-add contraction: faithful mode reproduces the original programs exactly
target_compile_options(znmno2_f PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffp-contract=off -std=f2008 -Wall -Wno-unused-dummy-argument>)
