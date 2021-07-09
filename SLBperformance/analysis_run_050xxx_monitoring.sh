#!/bin/bash

run="run_050001"
run_file="run_050001_25062021_19h25min_Ascii.dat"
data_folder="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB2021/run_050XXX/run_050001_25062021_19h25min_Ascii"
output="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis/converter_SLB/convertedfiles/"${run}"/"
#mkdir $output

initial_folder=$PWD

for i in {0..2339}
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
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"retriggers\"\) &
    root -l -q DummyDisplay.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\"\)
done

cd results_monitoring
source hadd.sh "Monitoring_summary" $run
source hadd.sh "HitMapsSimpleTracks" $run
cd -

cd ../results_retriggers
source hadd.sh $run

