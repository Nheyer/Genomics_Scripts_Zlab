cmake_minimum_required(VERSION 3.14)

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
