cmake_minimum_required(VERSION 3.15)

project(Genomics_Scripts_Zlab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# argparse is header only, vendored as a submodule so we don't rely on a system install
set(ARGPARSE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External_tools/argparse/include)
if(NOT EXISTS ${ARGPARSE_INCLUDE_DIR}/argparse/argparse.hpp)
    message(FATAL_ERROR "argparse submodule is missing, run: git submodule update --init --recursive")
endif()

add_executable(MSAtoCONSENSUS Scripts/MSA_fasta_to_Consensus.cpp)
target_include_directories(MSAtoCONSENSUS PRIVATE ${ARGPARSE_INCLUDE_DIR})

# We build in source, so `clean` also takes out what configuring generated.
# The binary itself is already handled by the built in clean; we deliberately
# don't list bin/ so that anything else you left in there survives.
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt;\
${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles;\
${CMAKE_CURRENT_SOURCE_DIR}/Makefile;\
${CMAKE_CURRENT_SOURCE_DIR}/cmake_install.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json")
