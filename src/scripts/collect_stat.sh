#!/bin/bash

function collect_stat
{
  local task_name="$1"
  local gep_bin="$2"
  local config_file="$3"
  local input_samples_file="$4"
  local output_samples_file="$5"
  local output_stat_file="$6"
  local timeout="$7"
  local runs_count="$8"
  local diffs_count="$9"
  local mse_success_margin="${10}"

  #echo "task_name ${task_name}"
  #echo "gep_bin ${gep_bin}"
  #echo "config_file ${config_file}"
  #echo "input_samples_file ${input_samples_file}"
  #echo "output_samples_file ${output_samples_file}"
  #echo "output_stat_file ${output_stat_file}"
  #echo "timeout ${timeout}"
  #echo "runs_count ${runs_count}"
  #echo "diffs_count ${diffs_count}"
  #echo "mse_success_margin ${mse_success_margin}"

  local work_dir=$(dirname ${config_file})
  local work_prefix="${work_dir}/collect_stat"
  local tmp_log_file="${work_prefix}_log"
  local tmp_time_file="${work_prefix}_time"
  local tmp_filtered_log_file="${work_prefix}_filtered_log"
  local tmp_output_samples_file="${work_prefix}_output_samples"
  local tmp_output_best_samples_file="${work_prefix}_best_samples"
  local all_mses=()

  local short_task_name=$(printf "%-35s" ${task_name})

  local best_mse=10000000
  local total_time=0
  local successful_count=0

  for run_idx in `seq 1 ${runs_count}`
  do
    local progress=$(( ((${run_idx} - 1) * 100) / ${runs_count} ))
    echo -ne "${short_task_name}: ${progress}%\r"

    local time_per_current_run=${timeout}

    local gep_command_line="${gep_bin} ${config_file} ${input_samples_file} ${tmp_output_samples_file} 0"

    if [ ${diffs_count} -eq 0 ]; then
      if [ ${timeout} -eq 0 ]; then
        /usr/bin/time -f "%e" -o ${tmp_time_file} ${gep_command_line} > ${tmp_log_file}
        time_per_current_run=$(cat ${tmp_time_file})
      else
        timeout --signal=SIGUSR1 ${timeout} ${gep_command_line} > ${tmp_log_file}
      fi
      total_time=$(echo "scale=2; ${total_time} + ${time_per_current_run}" | bc)
    else
      local tmp_input_samples_file=${input_samples_file}
      for diff_run_idx in `seq 1 10`
      do
        local prev_tmp_input_samples_file=${tmp_input_samples_file}
      done
    fi

    cat ${tmp_log_file} | grep "^| " | awk '{ print $2 " " $4 }' > ${tmp_filtered_log_file}
    local mse=`tail --lines=1 ${tmp_filtered_log_file} | awk '{ print $2 }'`

    all_mses[${run_idx}]=${mse}

    local better=$(echo ${mse} ${best_mse} | awk '{if ($1 <= $2) print "1"; else print "0"}')
    if [ ${better} -eq 1 ]; then
      best_mse=${mse}
      cp ${tmp_output_samples_file} ${tmp_output_best_samples_file}
    fi

    local successful=$(echo ${mse} ${mse_success_margin} | awk '{if ($1 <= $2) print "1"; else print "0"}')
    if [ ${successful} -eq 1 ]; then
      successful_count=$(( ${successful_count} + 1 ))
    fi
  done

  local avg_time_per_run=$(echo "scale=3; $total_time / ${runs_count}" | bc)
  local all_mses_line=${all_mses[*]}
  local success_rate_percentage=$(echo "scale=0; 100 * ${successful_count} / ${runs_count}" | bc)

  echo -ne "                                                     \r"
  printf "%-35s: best %4.3f; success %3d%% avg_time %4.1f\n" ${task_name} ${best_mse} ${success_rate_percentage} ${avg_time_per_run}

  echo "hostname `hostname`" > ${output_stat_file}
  echo "datetime `date`" >> ${output_stat_file}
  echo "runs_count ${runs_count}" >> ${output_stat_file}
  echo "timeout ${timeout}" >> ${output_stat_file}
  echo "avg_time ${avg_time_per_run}" >> ${output_stat_file}
  echo "mses ${all_mses_line}"  >> ${output_stat_file}
  echo "best_mse ${best_mse}"  >> ${output_stat_file}
  echo "mse_success_margin ${mse_success_margin}"  >> ${output_stat_file}
  echo "success_rate ${success_rate_percentage}"  >> ${output_stat_file}
  cat ${config_file} >> ${output_stat_file}

  cp ${tmp_output_best_samples_file} ${output_samples_file}
  rm -rf ${work_prefix}_*
}