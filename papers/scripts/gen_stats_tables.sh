#!/bin/bash

STATS_DIR=../stats
OUT_DIR=./build

function gen_plot
{
  for plot_name in "$@"; do
  done
}

gen_plot "pop_size_10" "pop_size_20" "pop_size_40" "pop_size_100" "pop_size_200" "pop_size_800"
