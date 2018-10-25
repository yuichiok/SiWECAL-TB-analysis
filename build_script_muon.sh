#!/bin/bash

conversion=$1
TDCrun=$2
path=$3
filename=$4

if [ -z "$4" ]; 
then 
    echo "Usage: source build_script.sh"
    echo " with 4 arguments needed: "
    echo "                    conversion =true or false "
    echo "                    TDCrun = true or false (if the TDC is enabled in the FEV13s or not)"
    echo "                    path (where the files are)"
    echo "                    file name without the dif_1_X_X.raw"
else
    echo " Running the building script "
    if [ $conversion == true ]
	then
	echo " Step 1 of 3: Converting the data files "
	if [ $TDCrun == false ]
	then
	    for dif in {"dif_1_1_1","dif_1_1_2","dif_1_1_3","dif_1_1_4","dif_1_1_5","dif_1_2_1","dif_1_2_2","dif_1_2_3","dif_1_2_4","dif_1_2_5"}
	    do
		root -l ConvertData.cc\(\"$path/$filename$dif.raw\"\)
	    done
	else 
	    for dif in {"dif_1_1_1","dif_1_1_2","dif_1_1_3","dif_1_1_4","dif_1_2_1","dif_1_2_2"}
            do
                root -l ConvertData.cc\(\"$path/$filename$dif.raw\"\)
            done
	    for dif in {"dif_1_1_5","dif_1_2_3","dif_1_2_4","dif_1_2_5"}
            do
                root -l ConvertDataTDC.cc\(\"$path/$filename$dif.raw\"\)
            done
	fi
	    
    else
        echo " Skipping step 1 of 3: conversion done before by hand "
    fi
fi
	


