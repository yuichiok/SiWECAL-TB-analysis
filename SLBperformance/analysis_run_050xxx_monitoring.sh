#!/bin/bash

run=$1
conversion=$2
shifter=$3
run_file="converted_run"
output=${PWD}"/../converter_SLB/convertedfiles/run_"${run}"/"

initial_folder=$PWD

if [ "$conversion" = true ]; then
    source conversion_run_050xxx.sh $run
fi

cd $initial_folder

j=0
for i in {0..999}
do
    
    if [ $i -gt 999 ]
    then
        i0=$i
    elif [ $i -gt 99 ]
    then
        i0=0$i
    elif [ $i -gt 9 ]
    then
        i0=00$i
    else
        i0=000$i
    fi

    if [ $j -gt 999 ]
    then
        i1=$j
    elif [ $j -gt 99 ]
    then
        i1=0$j
    elif [ $j -gt 9 ]
    then
        i1=00$j
    else
        i1=000$j
    fi
    FILE0=${output}"converted_run_"${run}".dat_"$i0".root"
    FILE1=${output}"converted_run_"${run}".dat_"$i1".root"
    echo $FILE0 $FILE1
    if [ -f "$FILE0" ]; then
        if [ -f "$FILE1" ]; then
	    cd $initial_folder
	    if [ "$shifter" = true ]; then
		root -l -q Monitoring.cc\(\"${output}/${run_file}_${run}.dat_$i0\",\"${run}_$j\",7,true\)
	    else
		root -l -q Monitoring.cc\(\"${output}/${run_file}_${run}.dat_$i0\",\"${run}_$j\",7,false\)
	    fi
	    cd -
	    j=$((i+1))
        else
            break
        fi
    fi
done

cd results_monitoring
source hadd.sh "Monitoring_summary" $run 
#source hadd.sh "HitMapsSimpleTracks" $run
cd -

