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
    echo " Step 2 of 3: Merging the root files"
    echo " #WARNING 1 !!"
    echo " MAKE ATTENTION TO THE lines 202 to 234 on the mergeRootFiles.cc file"
    echo " the factor x2 is applied to the slabs with 2.5MHz. For the others, we keep the default value of the bcid. Between the slabs at 5MHz and 2.5MHz, there is an offset. "
    ./mergeRootFiles.py $path/$filename ../built_files/${filename}_merge.root
    
    echo " Step 3 of 3: Building the events."
    echo " ----------------------------------------"
    echo " #WARNING 2 !!"
    echo " MAKE ATTENTION TO THE SLAB MAP defined in help_tools.py because"
    echo " it is using different files for mip and pedestal correction "
    echo " since for some of the slabs, this is not calculated."
    echo " ----------------------------------------"
    ./build_events.py ../built_files/${filename}_merge.root
fi
	


