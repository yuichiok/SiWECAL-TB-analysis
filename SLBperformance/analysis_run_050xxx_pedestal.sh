#!/bin/bash

#run
run="run_050010"
run_file="converted.dat"
#ascii folder
#rootfiles folder
output=${PWD}"/../converter_SLB/convertedfiles/run_050010_07172021_13h52min_Ascii/"

initial_folder=$PWD

for i in {1031..1774}
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
    #analysis
    cd $initial_folder
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"pedestal\"\) 
    cd ${initial_folder}/results_proto/
    if [ $i -gt 1 ]; then
	hadd Pedestal_new.root Pedestal_15_layers_${run}_${j}.root Pedestal_15_layers_${run}.root
	mv Pedestal_new.root Pedestal_15_layers_${run}.root
    else 
	mv Pedestal_15_layers_${run}_${j}.root Pedestal_15_layers_${run}.root
    fi
    tar czvf Pedestal_15_layers_${run}_${j}.root.tar.gz Pedestal_15_layers_${run}_${j}.root
    rm Pedestal_15_layers_${run}_${j}.root
    cd -

done

mv results_proto/Pedestal_15_layers_${run}.root results_pedestal/Pedestal_15_layers_${run}.root

cd results_pedestal
root -l -q PedestalFile.C\(\"Pedestal_${run}\"\)
cd ..

cp results_pedestal/Pedestal_${run}.txt ../pedestals/Pedestal_15_layers_${run}.txt

source analysis_run_050xxx_mips.sh
