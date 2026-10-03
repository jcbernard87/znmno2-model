add_executable(znmno2_cpp znmno2.cpp)
# No fused multiply-add contraction: faithful mode reproduces the original programs exactly
target_compile_options(znmno2_cpp PRIVATE $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:-ffp-contract=off -Wall -Wextra>)
