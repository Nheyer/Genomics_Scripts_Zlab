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

# --- tests -------------------------------------------------------------------
# CUnit is vendored as a submodule too. We only want the library, not CUnit's
# own test suite and examples.
set(CUNIT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External_tools/cunit)
if(EXISTS ${CUNIT_DIR}/CMakeLists.txt)
    set(CUNIT_DISABLE_TESTS TRUE CACHE BOOL "" FORCE)
    set(CUNIT_DISABLE_EXAMPLES TRUE CACHE BOOL "" FORCE)
    add_subdirectory(${CUNIT_DIR} EXCLUDE_FROM_ALL)

    enable_testing()
    add_executable(test_consensus Tests/test_consensus.cpp)
    # Tests/ includes the program source directly, see the comment at the top
    # of test_consensus.cpp for why.
    target_include_directories(test_consensus PRIVATE
        ${ARGPARSE_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/Scripts)
    target_link_libraries(test_consensus PRIVATE cunit)
    add_test(NAME consensus COMMAND test_consensus)
else()
    message(STATUS "cunit submodule missing, skipping tests. "
                   "Run: git submodule update --init --recursive")
endif()

# We build in source, so `clean` also takes out what configuring generated.
# The binary itself is already handled by the built in clean; we deliberately
# don't list bin/ so that anything else you left in there survives.
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt;\
${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles;\
${CMAKE_CURRENT_SOURCE_DIR}/Makefile;\
${CMAKE_CURRENT_SOURCE_DIR}/cmake_install.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json;\
${CMAKE_CURRENT_SOURCE_DIR}/CTestTestfile.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/Testing;\
${CMAKE_CURRENT_SOURCE_DIR}/CPackConfig.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/CPackSourceConfig.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/_CPack_Packages;\
${CUNIT_DIR}/CMakeFiles;\
${CUNIT_DIR}/Makefile;\
${CUNIT_DIR}/cmake_install.cmake;\
${CUNIT_DIR}/CTestTestfile.cmake;\
${CUNIT_DIR}/CUnit/CMakeFiles;\
${CUNIT_DIR}/CUnit/Makefile;\
${CUNIT_DIR}/CUnit/cmake_install.cmake;\
${CUNIT_DIR}/CUnit/CTestTestfile.cmake")
