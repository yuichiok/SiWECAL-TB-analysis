#!/bin/bash

run=$1
conversion=$2
run_file="converted"
output=${PWD}"/../converter_SLB/convertedfiles/"${run}"/"

PEDESTAL=1
MIP=1
MONITORING=0

initial_folder=$PWD

if [ $conversion -gt 0 ]; then
    source conversion_run_050xxx.sh ${run}
fi

if [ $MONITORING -gt 0 ]; then
    source analysis_run_050xxx_monitoring.sh ${run} 0
fi


cd $initial_folder

j=0

if [ $PEDESTAL -gt 0 ]; then
    source analysis_run_050xxx_pedestal.sh ${run}
fi

if [ $MIP -gt 0 ]; then
    source analysis_run_050xxx_mips.sh ${run}
fi
