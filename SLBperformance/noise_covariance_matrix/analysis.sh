initial=${PWD}
cd ../../converter_SLB/convertedfiles/
for run in "fromMIPScan_LowEnergyElectrons"
do

    cd ${initial}/analysis
    for sca in {1..15}
    do
	# run="03102022_pedestal_13slabs"
	root -l -q NoiseStudy.C\(\"${run}\",${sca}\)
    done
    root -l -q SummaryPlots.C\(\"$run\"\)
    cd ${initial}
done



	
