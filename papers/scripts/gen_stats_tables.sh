#!/bin/bash

STATS_DIR=../stats
OUT_DIR=./build

function gen_plot
{
  local ref_name="$1"
  local out_file="${OUT_DIR}/${ref_name}.tex"
  local tbl_label="tbl:${ref_name}"
  local first_file="${STATS_DIR}/$2_sin"

  local caption=$(grep "test_caption" ${first_file} | cut -d' ' -f2- | sed 's/^ *//')

  echo "
\begin{table}[!h]
  \caption{${caption}}
  \label{${tbl_label}}
  \begin{center}
    \begin{tabular}{|l|c|c|c|c|c|c|}
      \hline
      \multirow{2}{*}{Эксперимент} & \multicolumn{2}{ c| }{\$sin\$} & \multicolumn{2}{ c| }{\$rosen\$} & \multicolumn{2}{ c| }{\$sigmas\$} \\\\
      \cline{2-7}" > ${out_file}

echo "      & \$e_{b}\$ & \$r_{f}\$ & \$e_{b}\$ & \$r_{f}\$ & \$e_{b}\$ & \$r_{f}\$ \\\\
      \hline" >> ${out_file}

  local args_to_skip=1
  for arg in "$@"; do
    if [ ! ${args_to_skip} -eq 0 ]; then
      args_to_skip=$(( ${args_to_skip} - 1))
      continue
    fi

    local stat_file_prefix="${STATS_DIR}/${arg}"

    local stat_file_sin="${stat_file_prefix}_sin"
    local stat_file_rosen="${stat_file_prefix}_rosen"
    local stat_file_ogr1="${stat_file_prefix}_ogr1"

    local best_mse_sin=$(grep "best_mse" ${stat_file_sin} | awk '{ print $2 }')
    local best_mse_rosen=$(grep "best_mse" ${stat_file_rosen} | awk '{ print $2 }')
    local best_mse_ogr1=$(grep "best_mse" ${stat_file_ogr1} | awk '{ print $2 }')

    local success_rate_sin=$(grep "success_rate" ${stat_file_sin} | awk '{ print $2 }')
    local success_rate_rosen=$(grep "success_rate" ${stat_file_rosen} | awk '{ print $2 }')
    local success_rate_ogr1=$(grep "success_rate" ${stat_file_ogr1} | awk '{ print $2 }')

    local experiment_caption=$(grep "experiment_caption" ${stat_file_sin} | cut -d' ' -f2- | sed 's/^ *//')

    echo "      ${experiment_caption} & ${best_mse_sin} & ${success_rate_sin} & ${best_mse_rosen} & ${success_rate_rosen} & ${best_mse_ogr1} & ${success_rate_ogr1} \\\\" >> ${out_file}
  done

  echo "      \hline
    \end{tabular}
  \end{center}
\end{table}
" >> ${out_file}
}

gen_plot "cmp_pop_sizes" "pop_size_10" "pop_size_20" "pop_size_40" "pop_size_100" "pop_size_200" "pop_size_800"

gen_plot "cmp_mutation_probabilites" \
         "mutations_0_ferreira"  "mutations_1_ferreira"  "mutations_2_ferreira"  "mutations_3_ferreira"  "mutations_5_ferreira"  "mutations_10_ferreira" \
         "mutations_0_prefix"  "mutations_1_prefix"  "mutations_2_prefix"  "mutations_3_prefix"  "mutations_5_prefix"  "mutations_10_prefix" \
         "mutations_0_overlapped"  "mutations_1_overlapped"  "mutations_2_overlapped"  "mutations_3_overlapped"  "mutations_5_overlapped"  "mutations_10_overlapped"

gen_plot "cmp_dynamic_constants" "dynamic_constants" "dyn_smooth_constants"

gen_plot "cmp_tree_depths" \
         "tree_depth_3_ferreira" "tree_depth_4_ferreira" "tree_depth_5_ferreira" "tree_depth_6_ferreira" "tree_depth_7_ferreira" \
         "tree_depth_3_prefix" "tree_depth_4_prefix" "tree_depth_5_prefix" "tree_depth_6_prefix" "tree_depth_7_prefix" \
         "tree_depth_3_overlapped" "tree_depth_4_overlapped" "tree_depth_5_overlapped" "tree_depth_6_overlapped" "tree_depth_7_overlapped"

gen_plot "cmp_probability_density_selection_and_tournament" \
         "sel_roulette_10_ferreira" "sel_roulette_50_ferreira" "sel_roulette_400_ferreira" \
         "sel_prob_dens_10_ferreira" "sel_prob_dens_50_ferreira" "sel_prob_dens_400_ferreira" \
         "sel_tournament_10_ferreira" "sel_tournament_50_ferreira" "sel_tournament_400_ferreira" \
         "sel_roulette_10_prefix" "sel_roulette_50_prefix" "sel_roulette_400_prefix" \
         "sel_prob_dens_10_prefix" "sel_prob_dens_50_prefix" "sel_prob_dens_400_prefix" \
         "sel_tournament_10_prefix" "sel_tournament_50_prefix" "sel_tournament_400_prefix" \
         "sel_roulette_10_overlapped" "sel_roulette_50_overlapped" "sel_roulette_400_overlapped" \
         "sel_prob_dens_10_overlapped" "sel_prob_dens_50_overlapped" "sel_prob_dens_400_overlapped" \
         "sel_tournament_10_overlapped" "sel_tournament_50_overlapped" "sel_tournament_400_overlapped"

gen_plot "cmp_mse_k" \
         "mse_k_001_ferreira" "mse_k_01_ferreira" "mse_k_1_ferreira" "mse_k_10_ferreira" "mse_k_100_ferreira" \
         "mse_k_001_prefix" "mse_k_01_prefix" "mse_k_1_prefix" "mse_k_10_prefix" "mse_k_100_prefix" \
         "mse_k_001_overlapped" "mse_k_01_overlapped" "mse_k_1_overlapped" "mse_k_10_overlapped" "mse_k_100_overlapped"

gen_plot "cmp_replace_worst" "no_replace_worst" "replace_worst"

gen_plot "cmp_fitnesses" \
         "full_mse_ferreira" "partial_mse_ferreira" "r_squared_ferreira" \
         "full_mse_prefix" "partial_mse_prefix" "r_squared_prefix" \
         "full_mse_overlapped" "partial_mse_overlapped" "r_squared_overlapped"

gen_plot "cmp_additional_population" "no_add_pop" "add_pop"
