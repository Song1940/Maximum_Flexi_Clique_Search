cmake_minimum_required(VERSION 3.15)
project(FlexiClique CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(FLEXI_NATIVE "Compile with -march=native, as for the experiments of the paper" ON)
option(FLEXI_VARIANTS "Build the variant binaries of the experiments (build/variants/)" ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(FLEXI_NATIVE)
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native")
endif()

set(FLEXI_SOURCES
  src/main.cpp src/graph.cpp src/fpa.cpp src/peel_split.cpp src/eba.cpp src/npa.cpp)

# FPA, NPA, and EBA.
add_executable(flexi ${FLEXI_SOURCES})
target_include_directories(flexi PRIVATE src)
target_compile_options(flexi PRIVATE -Wall -Wextra)

# Variants of the experiments: the same sources with the macros documented in
# src/eba.cpp. Each takes the options of flexi.
#   name             macros                                            experiment
set(FLEXI_VARIANT_LIST
  "no_rule1         NO_RULE1"                                        # EQ4 (EBA-Rule 1)
  "no_rule2         NO_RULE2"                                        # EQ4 (EBA-Rule 2)
  "no_rule3         NO_RULE3"                                        # EQ4 (EBA-Rule 3)
  "no_incumbent     NO_INCUMBENT"                                    # EQ4 (EBA-FPA)
  "no_pivot         NO_PIVOT NO_CANDIDATE_SATURATION"                # EQ5 (EBA-Pivot)
  "no_degeneracy    NO_DEGENERACY_ORDER"                             # EQ5 (EBA-Degeneracy)
  "binary           BINARY_BRANCHING NO_CANDIDATE_SATURATION"        # EQ5 (Binary)
  "flat_symse       FLAT_SYMSE NO_CANDIDATE_SATURATION"              # EQ5 (flat Sym-SE)
  "root_stats       ROOT_STATS"                                      # EQ2 (nodes left after the root pruning)
  "width            TRACE_WIDTH"                                     # technical report: branching width
  "width_no_pivot   TRACE_WIDTH NO_PIVOT NO_CANDIDATE_SATURATION"
)
if(FLEXI_VARIANTS)
  foreach(entry ${FLEXI_VARIANT_LIST})
    string(REGEX REPLACE " +" ";" fields "${entry}")
    list(POP_FRONT fields name)
    add_executable(flexi_${name} ${FLEXI_SOURCES})
    target_include_directories(flexi_${name} PRIVATE src)
    target_compile_options(flexi_${name} PRIVATE -Wall -Wextra)
    target_compile_definitions(flexi_${name} PRIVATE ${fields} FLEXI_VARIANT="${name}")
    set_target_properties(flexi_${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/variants)
  endforeach()
endif()

# Baseline of the case study (EQ7): maximum clique, k-plex, or quasi-clique
# containing a query node.
add_executable(qcs baselines/qcs/qcs.cpp src/graph.cpp)
target_include_directories(qcs PRIVATE src)
target_compile_options(qcs PRIVATE -Wall -Wextra)
