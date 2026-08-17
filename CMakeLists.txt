cmake_minimum_required(VERSION 3.14)

project(Genomics_Scripts_Zlab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable(MSAtoCONSENSUS Scripts/MSA_fasta_to_Consensus.cpp)
