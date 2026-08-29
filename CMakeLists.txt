cmake_minimum_required(VERSION 3.15)

project(Genomics_Scripts_Zlab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# argparse is header only, vendored as a submodule so we don't rely on a system install
set(ARGPARSE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External_tools/argparse/include)
if(NOT EXISTS ${ARGPARSE_INCLUDE_DIR}/argparse/argparse.hpp)
    message(FATAL_ERROR "argparse submodule is missing, run: git submodule update --init --recursive")
endif()

# Each tool lives in its own directory under CppSrc/, which holds source only.
# Data (the enzyme table, the test fixtures) lives under Data/.
set(MSA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CppSrc/MSA-to-consensus)
set(DIGEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CppSrc/Enzyme-digest)
set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Data)

# Target names match the directory names, and the binaries in bin/ get the
# same names again.
add_executable(MSA-to-consensus ${MSA_DIR}/MSA_fasta_to_Consensus.cpp)
target_include_directories(MSA-to-consensus PRIVATE ${ARGPARSE_INCLUDE_DIR})

# --- enzyme table ------------------------------------------------------------
# Data/restriction_enzymes.csv is the single source of truth for the enzyme
# data, so it is embedded verbatim rather than re-encoded as C++ literals: two
# copies would drift. The generated header goes in the build tree, never the
# source tree, and the same parser reads it that --enzyme-db uses on a file.
# Re-run cmake after editing the CSV to pick the change up.
set(GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(READ ${DATA_DIR}/restriction_enzymes.csv ENZYME_CSV_CONTENTS)
configure_file(
    ${DIGEST_DIR}/enzyme_csv_embedded.hpp.in
    ${GENERATED_INCLUDE_DIR}/enzyme_csv_embedded.hpp
    @ONLY)

add_executable(Enzyme-digest ${DIGEST_DIR}/Restriction_Enzyme_Digest.cpp)
target_include_directories(Enzyme-digest PRIVATE
    ${ARGPARSE_INCLUDE_DIR}
    ${GENERATED_INCLUDE_DIR}
    ${DIGEST_DIR})

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
        ${MSA_DIR})
    target_link_libraries(test_consensus PRIVATE cunit)
    add_test(NAME consensus COMMAND test_consensus)

    add_executable(test_digest Tests/test_digest.cpp)
    target_include_directories(test_digest PRIVATE
        ${ARGPARSE_INCLUDE_DIR}
        ${GENERATED_INCLUDE_DIR}
        ${DIGEST_DIR})
    target_link_libraries(test_digest PRIVATE cunit)
    add_test(NAME digest COMMAND test_digest)
else()
    message(STATUS "cunit submodule missing, skipping tests. "
                   "Run: git submodule update --init --recursive")
endif()

# --- install -----------------------------------------------------------------
# Only the tools go out; test_consensus is for us, and cunit is a build
# dependency that EXCLUDE_FROM_ALL above already keeps out of the install.
include(GNUInstallDirs)
install(TARGETS MSA-to-consensus Enzyme-digest
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# LICENSE.restriction-digest is the MIT notice for the code Enzyme-digest was
# ported from; it has to travel with the binary we ship.
install(FILES README.md LICENSE LICENSE.restriction-digest
    DESTINATION ${CMAKE_INSTALL_DOCDIR})

# We build in source, so `clean` also takes out what configuring generated.
# The binary itself is already handled by the built in clean; we deliberately
# don't list bin/ so that anything else you left in there survives.
set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeCache.txt;\
${CMAKE_CURRENT_SOURCE_DIR}/CMakeFiles;\
${GENERATED_INCLUDE_DIR};\
${CMAKE_CURRENT_SOURCE_DIR}/Makefile;\
${CMAKE_CURRENT_SOURCE_DIR}/cmake_install.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json;\
${CMAKE_CURRENT_SOURCE_DIR}/CTestTestfile.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/Testing;\
${CMAKE_CURRENT_SOURCE_DIR}/CPackConfig.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/CPackSourceConfig.cmake;\
${CMAKE_CURRENT_SOURCE_DIR}/_CPack_Packages;\
${CMAKE_CURRENT_SOURCE_DIR}/install_manifest.txt;\
${CUNIT_DIR}/CMakeFiles;\
${CUNIT_DIR}/Makefile;\
${CUNIT_DIR}/cmake_install.cmake;\
${CUNIT_DIR}/CTestTestfile.cmake;\
${CUNIT_DIR}/CUnit/CMakeFiles;\
${CUNIT_DIR}/CUnit/Makefile;\
${CUNIT_DIR}/CUnit/cmake_install.cmake;\
${CUNIT_DIR}/CUnit/CTestTestfile.cmake")
