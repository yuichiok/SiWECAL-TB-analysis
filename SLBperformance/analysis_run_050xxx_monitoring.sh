#!/bin/bash

run="run_050004"
run_file="converted.dat"
output=${PWD}"/../converter_SLB/convertedfiles/run_050004_09072021_21h05min_Ascii/"


initial_folder=$PWD

for i in {1..1803}
do
    j="000"$i
    if [ $i -gt 9 ]; then
        j="00"$i
    fi
    
    if [ $i -gt 99 ]; then
        j="0"$i
    fi

    if [ $i -gt 999 ]; then
        j=$i
    fi
    #conversion
    #analysis
    cd $initial_folder
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"retriggers\"\)  &
    root -l -q DummyDisplay.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\"\)
done

cd results_monitoring
source hadd.sh "Monitoring_summary" $run &
source hadd.sh "HitMapsSimpleTracks" $run
cd -

cd ../results_retriggers
source hadd.sh $run

