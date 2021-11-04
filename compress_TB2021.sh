#!/bin/bash

run=$1
folder="/mnt/win1/Run_data/run_"$run
cd $folder

tar czvf "run_"${run}".dat".tar.gz "run_"${run}".dat"

	    
j=1
for i in {0..9999}
do
    FILE_start="run_"${run}".dat"
    if test -f "${FILE_start}"; then
	FILE_new=${FILE_start}"_0000"
	mv ${FILE_start} ${FILE_new}
	tar czvf ${FILE_new}.tar.gz ${FILE_new}
    fi
    j=$((i+1))

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

    echo $i0 $i1

    FILE0="run_"${run}".dat_"$i0
    FILE1="run_"${run}".dat_"$i1
    if [ -f "$FILE0" ]; then
	if [ -f "$FILE1" ]; then
	    tar czvf ${FILE0}.tar.gz $FILE0
	    rm $FILE0
	else
	    break
	fi
    fi
done

cd -
