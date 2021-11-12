#!/bin/bash

run=$1
folder="/mnt/win1/Run_data/"$run
cd $folder

	    
j=1
for i in {0..9999}
do
    FILE_start=${run}".dat"
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

    #echo $i0 

    FILE0=${run}".dat_"$i0
    if [ -f "$FILE0" ]; then
	tar czvf ${FILE0}.tar.gz $FILE0
	rm $FILE0
    fi

done

cd - 
#mv $folder /mnt/win1/Run_data/compressed/.

