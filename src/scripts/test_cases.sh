#!/bin/bash

. collect_stat.sh

function create_gep_config_file
{
  local out_cfg_file="$1"

  echo -n > ${out_cfg_file}
  echo "test_caption       ${TEST_CAPTION}"  >> ${out_cfg_file}
  echo "experiment_caption ${EXPERIMENT_CAPTION}" >> ${out_cfg_file}
  echo "population_size    ${POP_SIZE}"      >> ${out_cfg_file}
  echo "tree_depth         ${TREE_DEPTH}"    >> ${out_cfg_file}
  echo "genes_count        ${GENES_COUNT}"   >> ${out_cfg_file}
  echo "generations_count  ${GENERATIONS_COUNT}"  >> ${out_cfg_file}
  echo "mutations_per_chromosome  ${MUTATIONS_PER_CHROMOSOME}" >> ${out_cfg_file}
  echo "ops_preset         ${OPS_PRESET}"    >> ${out_cfg_file}
  echo "coding_type        ${CODING_TYPE}"   >> ${out_cfg_file}
  echo "fitness_type       ${FITNESS_TYPE}"  >> ${out_cfg_file}
  echo "mse_coefficient    ${MSE_FITNESS_K}" >> ${out_cfg_file}

  if [ $USE_TOURNAMENT -eq 1 ]; then
    echo "use_tournament 1" >> ${out_cfg_file}
  fi

  if [ $USE_REPLACE_WORST -eq 1 ]; then
    echo "use_replace_worst 1" >> ${out_cfg_file}
  fi

  if [ $USE_INCREMENTAL_EVOLUTION -eq 1 ]; then
    echo "use_incremental_evolution 1" >> ${out_cfg_file}
  fi

  if [ $USE_SELECT_PROB_DENSITY -eq 1 ]; then
    echo "use_selection_probability_density 1" >> ${out_cfg_file}
  fi

  if [ $USE_ADDITIONAL_POPULATION -eq 1 ]; then
    echo "use_additional_population 1" >> ${out_cfg_file}
  fi

  if [ $USE_DYNAMIC_CONSTANTS -eq 1 ]; then
    echo "use_dynamic_constants 1" >> ${out_cfg_file}
  fi

  if [ $ENABLE_GRAPHICS -eq 1 ]; then
    echo "enable_graphics 1" >> ${out_cfg_file}
  fi
}

function set_origin_params
{
  POP_SIZE=100
  TREE_DEPTH=5
  GENES_COUNT=1
  GENERATIONS_COUNT=100000
  MUTATIONS_PER_CHROMOSOME=2

  OPS_PRESET="short"
  CODING_TYPE="ferreira"
  FITNESS_TYPE="mse"
  MSE_FITNESS_K="1.0"

  USE_TOURNAMENT=0
  USE_REPLACE_WORST=0
  USE_INCREMENTAL_EVOLUTION=0
  USE_SELECT_PROB_DENSITY=0
  USE_ADDITIONAL_POPULATION=0
  USE_DYNAMIC_CONSTANTS=0

  DIFFERENTIALS_COUNT=0

  MSE_SUCCESS_MARGIN="0.01"

  NO_TIME_LIMIT=0
  TIME_LIMIT_MILTIPLIER=1
}

function set_default_params
{
  set_origin_params
  MSE_FITNESS_K="10"
}

function set_task_1_params
{
  TEST_BASENAME="sin"
  TEST_FILE="${TEST_FILES_DIR}/test_sin"
  TIMEOUT=7
  MSE_SUCCESS_MARGIN=".08"
  OPS_PRESET="short"

  if [ ${NO_TIME_LIMIT} -eq 1 ]; then
    GENERATIONS_COUNT=4000
  fi
}

function set_task_2_params
{
  TEST_BASENAME="rosen"
  TEST_FILE="${TEST_FILES_DIR}/test_Rosenbrock"
  TIMEOUT=10
  MSE_SUCCESS_MARGIN=".05"
  OPS_PRESET="short"

  if [ ${NO_TIME_LIMIT} -eq 1 ]; then
    GENERATIONS_COUNT=1000
  fi
}

function set_task_3_params
{
  TEST_BASENAME="ogr1"
  TEST_FILE="${TEST_FILES_DIR}/test_OGR1"
  TIMEOUT=10
  MSE_SUCCESS_MARGIN=".06"
  OPS_PRESET="full"

  if [ ${NO_TIME_LIMIT} -eq 1 ]; then
    GENERATIONS_COUNT=1000
  fi
}

function run_task
{
  local test_full_name="$1"
  local config_file="$2"

  # Set short timeout for quick testing
  #TIMEOUT=1

  # Set no time limit - process will finish when all set generations are processed
  if [ ${NO_TIME_LIMIT} -eq 1 ]; then
    TIMEOUT=0
  fi

  TIMEOUT=$((${TIMEOUT} * ${TIME_LIMIT_MILTIPLIER}))

  create_gep_config_file ${PATH_TO_GENERATED_CONFIG_FILE}

  collect_stat "${test_full_name}" "${GEP_TEST_BIN}" "${config_file}" "${TEST_FILE}" "${WORK_DIR}/needless_output_samples" "${STAT_OUTPUT_DIR}/${test_full_name}" ${TIMEOUT} ${COLLECT_RUNS_COUNT} ${DIFFERENTIALS_COUNT} "${MSE_SUCCESS_MARGIN}"
}

function run_all_tasks
{
  EXPERIMENT_CAPTION="$1"
  local test_name="$2"

  set_task_1_params
  run_task "${test_name}_${TEST_BASENAME}" ${PATH_TO_GENERATED_CONFIG_FILE}
  
  if [ ! ${JUST_FIRST_TASK} -eq 1 ]; then
    set_task_2_params
    run_task "${test_name}_${TEST_BASENAME}" ${PATH_TO_GENERATED_CONFIG_FILE}

    set_task_3_params
    run_task "${test_name}_${TEST_BASENAME}" ${PATH_TO_GENERATED_CONFIG_FILE}
  fi
}

function run_all_tasks_all_codings
{
  local experiment_caption_prefix="$1"
  local test_name="$2"

  CODING_TYPE="ferreira"
  run_all_tasks "${experiment_caption_prefix} (ферр.)" "${test_name}_ferreira"

  CODING_TYPE="prefix"
  run_all_tasks "${experiment_caption_prefix} (префикс.)" "${test_name}_prefix"

  CODING_TYPE="overlapped"
  run_all_tasks "${experiment_caption_prefix} (налож.)" "${test_name}_overlapped"
}



function gep_test_origin
{
  TEST_CAPTION="Тестовый запуск"
  set_origin_params

  NO_TIME_LIMIT=1
  GENERATIONS_COUNT=100
  run_all_tasks "Исходный алгоритм" "origin"
}

function gep_test_pop_size
{
  TEST_CAPTION="Эффективность алгоритма при различных размерах популяции"
  set_default_params

  POP_SIZE=10;  run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=20;  run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=40;  run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=60;  run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=80;  run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=100; run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=120; run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=140; run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=200; run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
  POP_SIZE=800; run_all_tasks "${POP_SIZE} особей" "pop_size_${POP_SIZE}"
}

function gep_test_mutations
{
  TEST_CAPTION="Эффективность алгоритма при различных вероятностях мутации"
  set_default_params

  MUTATIONS_PER_CHROMOSOME=0;  run_all_tasks_all_codings "Без мутаций"                         "mutations_${MUTATIONS_PER_CHROMOSOME}"
  MUTATIONS_PER_CHROMOSOME=1;  run_all_tasks_all_codings "${MUTATIONS_PER_CHROMOSOME} мутаций" "mutations_${MUTATIONS_PER_CHROMOSOME}"
  MUTATIONS_PER_CHROMOSOME=2;  run_all_tasks_all_codings "${MUTATIONS_PER_CHROMOSOME} мутаций" "mutations_${MUTATIONS_PER_CHROMOSOME}"
  MUTATIONS_PER_CHROMOSOME=3;  run_all_tasks_all_codings "${MUTATIONS_PER_CHROMOSOME} мутаций" "mutations_${MUTATIONS_PER_CHROMOSOME}"
  MUTATIONS_PER_CHROMOSOME=5;  run_all_tasks_all_codings "${MUTATIONS_PER_CHROMOSOME} мутаций" "mutations_${MUTATIONS_PER_CHROMOSOME}"
}

function gep_test_dynamic_constants
{
  TEST_CAPTION="Эффективность алгоритма при плавной и динамической мутации констант"
  set_default_params

  USE_DYNAMIC_CONSTANTS=0; run_all_tasks "Плавно-динамическая" "dyn_smooth_constants"
  USE_DYNAMIC_CONSTANTS=1; run_all_tasks "Динамическая"        "dynamic_constants"
}

function gep_test_tree_depth
{
  TEST_CAPTION="Эффективность алгоритма при разной глубине синтаксического дерева"
  set_default_params

  TREE_DEPTH=3; run_all_tasks_all_codings "Глубина ${TREE_DEPTH}" "tree_depth_${TREE_DEPTH}"
  TREE_DEPTH=4; run_all_tasks_all_codings "Глубина ${TREE_DEPTH}" "tree_depth_${TREE_DEPTH}"
  TREE_DEPTH=5; run_all_tasks_all_codings "Глубина ${TREE_DEPTH}" "tree_depth_${TREE_DEPTH}"
  TREE_DEPTH=6; run_all_tasks_all_codings "Глубина ${TREE_DEPTH}" "tree_depth_${TREE_DEPTH}"
  TREE_DEPTH=7; run_all_tasks_all_codings "Глубина ${TREE_DEPTH}" "tree_depth_${TREE_DEPTH}"
}

function gep_test_selections
{
  TEST_CAPTION="Эффективность алгоритма с разными операторах отбора и размере популяции"
  set_default_params

  POP_SIZE=10;  run_all_tasks_all_codings "Рулетка ${POP_SIZE} особей" "sel_roulette_${POP_SIZE}"
  POP_SIZE=50;  run_all_tasks_all_codings "Рулетка ${POP_SIZE} особей" "sel_roulette_${POP_SIZE}"
  POP_SIZE=100  run_all_tasks_all_codings "Рулетка ${POP_SIZE} особей" "sel_roulette_${POP_SIZE}"
  POP_SIZE=150; run_all_tasks_all_codings "Рулетка ${POP_SIZE} особей" "sel_roulette_${POP_SIZE}"
  POP_SIZE=400; run_all_tasks_all_codings "Рулетка ${POP_SIZE} особей" "sel_roulette_${POP_SIZE}"

  USE_SELECT_PROB_DENSITY=1
  POP_SIZE=10;  run_all_tasks_all_codings "Отбор плотн. (${POP_SIZE})" "sel_prob_dens_${POP_SIZE}"
  POP_SIZE=50;  run_all_tasks_all_codings "Отбор плотн. (${POP_SIZE})" "sel_prob_dens_${POP_SIZE}"
  POP_SIZE=100; run_all_tasks_all_codings "Отбор плотн. (${POP_SIZE})" "sel_prob_dens_${POP_SIZE}"
  POP_SIZE=150; run_all_tasks_all_codings "Отбор плотн. (${POP_SIZE})" "sel_prob_dens_${POP_SIZE}"
  POP_SIZE=400; run_all_tasks_all_codings "Отбор плотню (${POP_SIZE})" "sel_prob_dens_${POP_SIZE}"
  USE_SELECT_PROB_DENSITY=0

  USE_TOURNAMENT=1
  POP_SIZE=10;  run_all_tasks_all_codings "Турнир ${POP_SIZE} особей" "sel_tournament_${POP_SIZE}"
  POP_SIZE=50;  run_all_tasks_all_codings "Турнир ${POP_SIZE} особей" "sel_tournament_${POP_SIZE}"
  POP_SIZE=100; run_all_tasks_all_codings "Турнир ${POP_SIZE} особей" "sel_tournament_${POP_SIZE}"
  POP_SIZE=150; run_all_tasks_all_codings "Турнир ${POP_SIZE} особей" "sel_tournament_${POP_SIZE}"
  POP_SIZE=400; run_all_tasks_all_codings "Турнир ${POP_SIZE} особей" "sel_tournament_${POP_SIZE}"
}

function gep_test_mse_k
{
  TEST_CAPTION="Эффективность алгоритма при различных коэффициентах давления \$K\$"
  set_default_params

  MSE_FITNESS_K="0.01"; run_all_tasks_all_codings "\$K=${MSE_FITNESS_K}\$" "mse_k_001"
  MSE_FITNESS_K="0.1";  run_all_tasks_all_codings "\$K=${MSE_FITNESS_K}\$" "mse_k_01"
  MSE_FITNESS_K="1";    run_all_tasks_all_codings "\$K=${MSE_FITNESS_K}\$" "mse_k_1"
  MSE_FITNESS_K="10";   run_all_tasks_all_codings "\$K=${MSE_FITNESS_K}\$" "mse_k_10"
  MSE_FITNESS_K="100";  run_all_tasks_all_codings "\$K=${MSE_FITNESS_K}\$" "mse_k_100"
}

function gep_test_replace_worst
{
  TEST_CAPTION="Эффективность алгоритма с заменой худших особей"
  set_default_params

  USE_REPLACE_WORST=0; run_all_tasks "Исходный алгоритм" "no_replace_worst"
  USE_REPLACE_WORST=1; run_all_tasks "Замена худших"     "replace_worst"
}

function gep_test_fitnesses
{
  TEST_CAPTION="Эффективность алгоритма с различными функциями фитнеса"
  set_default_params

  FITNESS_TYPE="mse";         run_all_tasks_all_codings "СКО"           "full_mse"
  FITNESS_TYPE="partial_mse"; run_all_tasks_all_codings "Частичная СКО" "partial_mse"
  FITNESS_TYPE="r_squared";   run_all_tasks_all_codings "\$R^2\$"       "r_squared"
}

function gep_test_incremental
{
  TEST_CAPTION="Эффективность алгоритма при традиционной и инкрементальной эволюции"
  set_default_params

  NO_TIME_LIMIT=1

  USE_INCREMENTAL_EVOLUTION=0; run_all_tasks "Традиц." "no_incremental"
  USE_INCREMENTAL_EVOLUTION=1; run_all_tasks "Инкрем." "incremental"
}

function gep_test_additional_population
{
  TEST_CAPTION="Эффективность алгоритма с дополнительной популяцией"
  set_default_params

  USE_ADDITIONAL_POPULATION=0; run_all_tasks "Исходный алгоритм" "no_add_pop"
  USE_ADDITIONAL_POPULATION=1; run_all_tasks "Доп. популяция"    "add_pop"
}

function gep_test_differential
{
  TEST_CAPTION="Эффективность алгоритма с разностным подходом"
  set_default_params

  TIME_LIMIT_MILTIPLIER=6; DIFFERENTIALS_COUNT=0; run_all_tasks "Исходный алгоритм" "no_differential"
  TIME_LIMIT_MILTIPLIER=3; DIFFERENTIALS_COUNT=1; run_all_tasks "Разность 1"        "differential_${DIFFERENTIALS_COUNT}"
  TIME_LIMIT_MILTIPLIER=2; DIFFERENTIALS_COUNT=2; run_all_tasks "Разность 2"        "differential_${DIFFERENTIALS_COUNT}"
  TIME_LIMIT_MILTIPLIER=1; DIFFERENTIALS_COUNT=3; run_all_tasks "Разность 3"        "differential_${DIFFERENTIALS_COUNT}"
}
