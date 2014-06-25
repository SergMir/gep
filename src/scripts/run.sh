#!/bin/bash

WORK_DIR=../build
GEP_TEST_BIN=${WORK_DIR}/gep_test
GEP_COMPARATOR_BIN=${WORK_DIR}/gep_comparator
STAT_OUTPUT_DIR=../../stats
TEST_FILES_DIR=../test_files_gen/build

PATH_TO_GENERATED_CONFIG_FILE=${WORK_DIR}/config

COLLECT_RUNS_COUNT=100
ENABLE_GRAPHICS=0

JUST_FIRST_TASK=0

cd scripts

. test_cases.sh

action="$1"

case "${action}" in
  "origin" )
    gep_test_origin
    ;;
  "all_tests" )
    gep_test_pop_size
    gep_test_mutations
    gep_test_dynamic_constants
    gep_test_tree_depth
    gep_test_selections
    gep_test_mse_k
    gep_test_replace_worst
    gep_test_fitnesses
    #gep_test_incremental
    gep_test_additional_population
    #gep_test_differential
    ;;
  "valgrind" )
    COLLECT_RUNS_COUNT=1
    JUST_FIRST_TASK=1
    #GEP_TEST_BIN="valgrind --leak-check=full --show-reachable=yes --track-origins=yes --leak-resolution=high --gen-suppressions=yes ${GEP_TEST_BIN}"
    GEP_TEST_BIN="valgrind --leak-check=full --show-reachable=yes --track-origins=yes --leak-resolution=high --suppressions=./valgrind.suppresions ${GEP_TEST_BIN}"
    gep_test_origin
    ;;
  *)
    ${action}
    ;;
esac