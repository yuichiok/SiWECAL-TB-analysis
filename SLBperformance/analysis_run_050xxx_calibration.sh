#!/bin/bash

run=$1
conversion=$2
run_file="converted"
output=${PWD}"/../converter_SLB/convertedfiles/"${run}"/"

PEDESTALMIP=1
MONITORING=0

initial_folder=$PWD

if [ $conversion -gt 0 ]; then
    source conversion_run_050xxx.sh ${run}
fi


if [ $MONITORING -gt 0 ]; then
    cd $initial_folder
    source analysis_run_050xxx_stats.sh ${run} 0
fi

if [ $PEDESTALMIP -gt 0 ]; then
    cd $initial_folder
    source analysis_run_050xxx_pedestalmip2.sh ${run}
fi
