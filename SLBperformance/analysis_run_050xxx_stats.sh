/bash

run=$1
run_file="converted"
output=${PWD}"/../converter_SLB/convertedfiles/"${run}"/"


initial_folder=$PWD

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
    FILE0=${output}"converted_"${run}".dat_"$i0".root"
    FILE1=${output}"converted_"${run}".dat_"$i1".root"
    echo $FILE0 $FILE1
    if [ -f "$FILE0" ]; then
        if [ -f "$FILE1" ]; then
	    cd $initial_folder
            if [ $((j%2)) -eq 0 ]
            then
		root -l -q Proto.cc\(\"${output}/${run_file}_${run}.dat_$i0\",\"${run}_$i0\",1\) &
            else
		root -l -q Proto.cc\(\"${output}/${run_file}_${run}.dat_$i0\",\"${run}_$i0\",1\)  
            fi
	    cd -
	    j=$((i+1))
        else
	    break
        fi
    fi
done

sleep 45
cd TBchecks/
j=$((j-1))
root -l -q statistics.cc\(\"HitMonitoring_${run}\",$j\)
rm results_proto/HitMonitoring_${run}*
cd -

