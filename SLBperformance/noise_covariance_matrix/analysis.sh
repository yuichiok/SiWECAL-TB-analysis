for sca in {0..14}
do
    # run="03102022_pedestal_13slabs"
    run="Run_ILC_03112022_pedestal"
    cd ../
    root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",0,${sca}\)
    root -l -q Noise.cc\(\"../converter_SLB/convertedfiles/${run}.root\",\"${run}\",1,${sca}\)
    cd -
done


	
