initial=${PWD}
cd ../../converter_SLB/convertedfiles/
for run in "Run_ILC_15052022_cosmic_it17"
do

    cd ${initial}/analysis
    for sca in {1..15}
    do
	# run="03102022_pedestal_13slabs"
	root -l -q NoiseStudy.C\(\"${run}\",${sca}\) &
    done
    root -l -q SummaryPlots.C\(\"$run\"\)
    cd ${initial}
done



	
