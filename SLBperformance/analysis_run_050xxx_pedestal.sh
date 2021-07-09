#!/bin/bash

#run
run="run_050001"
run_file="run_050001_25062021_19h25min_Ascii.dat"
#ascii folder
data_folder="/lustre/ific.uv.es/prj/ific/flc/SiWECAL/TB2021/run_050XXX/run_050001_25062021_19h25min_Ascii"
#rootfiles folder
output="/mnt/HardDrive/cernbox_hd/SiWECAL/TB2021/SiWECAL-TB-analysis/converter_SLB/convertedfiles/"${run}"/"
mkdir $output

initial_folder=$PWD

for i in {1..2339}
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
    #cd ${output}/
    #tar xzvf ${data_folder}/${run_file}_$j.tar.gz
    #cd $initial_folder/../converter_SLB
    #root -l -q ConvertDataSL.cc\(\"${output}/${run_file}_$j\",false\)
    #rm ${data_folder}/${run_file}_$j
    
    #analysis
    cd $initial_folder
    root -l -q Proto.cc\(\"${output}/${run_file}_$j\",\"${run}_$j\",\"pedestal\"\) 
    cd ${initial_folder}/results_proto/
    hadd Pedestal_new.root Pedestal_15_layers_${run}_*.root Pedestal_15_layers_${run}.root
    mv Pedestal_new.root Pedestal_15_layers_${run}.root
    rm Pedestal_15_layers_${run}_*.root
    cd -

done

cd results_proto
for k in {0..2}
do
    for j in {0..9}
    do
	for i in {0..9}
	do
	    tar xzvf Pedestal_15_layers_run_${run}_${k}${j}0${i}.root.tar.gz
	    ##mv lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis/SLBperformance/results_proto/Pedestal_15_layers_${run}_${k}${j}0${i}.root .
	    tar czvf Pedestal_15_layers_${run}_${k}${j}0${i}.root &
	done
	for i in {10..99}
	do
	    tar xzvf Pedestal_15_layers_run_${run}_${k}${j}${i}.root.tar.gz
	    ##mv lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis/SLBperformance/results_proto/Pedestal_15_layers_${run}_${k}${j}${i}.root .
	    tar czvf Pedestal_15_layers_${run}_${k}${j}${i}.root.tar.gz Pedestal_15_layers_${run}_${k}${j}${i}.root &
	done
	
	hadd ../results_pedestal/Pedestal_${run}_${k}_${j}.root Pedestal_15_layers_${run}_${k}${j}*.root
	rm Pedestal_15_layers_${run}_${k}${j}*.root
	##	rm lhome/ific/a/airqui/SiWECAL/SiWECAL-TB-analysis/SLBperformance/results_proto/*_${k}${j}*.root
    done
done
cd ../

hadd results_pedestal/Pedestal_${run}.root hadd results_pedestal/Pedestal_${run}_*.root

cd results_pedestal
root -l -q PedestalFile.C\(\"Pedestal_${run}\"\)
cd ..
