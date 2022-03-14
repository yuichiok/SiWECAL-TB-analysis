#!/bin/bash

for folder in "run_050XXX/"
do
    cd $folder
    for f in run_*; do
	cd ${f}
	for file in *; 
	do
	    if [[ "$file" == *"tar.gz" ]];then
		printf '%s\n' "$file"
	    else 
		echo "Compress File ->" $file
		tar czvf ${file}.tar.gz ${file}
		rm ${file}
	    fi
	done
	cd ../
    done
    cd ../
done
