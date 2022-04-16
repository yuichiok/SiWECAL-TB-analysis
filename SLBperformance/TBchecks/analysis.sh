#sleep 20s
for mode in 2 
do
    #root -l -q analysis.cc\(\"MIPscan_1.2pF_0.8pF\",\"high\",false,true,$mode,\"run_050571_injection_merged_LowEnergyElectrons\"\)  &
    #root -l -q analysis.cc\(\"MIPscan_1.2pF_0.8pF\",\"low\",false,true,$mode,\"run_050571_injection_merged_LowEnergyElectrons\"\) &
    root -l -q analysis.cc\(\"MIPscan_6pF\",\"high\",false,true,$mode,\"run_050575_injection_merged_ILCEnergyElectrons\"\)
done
#root -l -q analysis.cc\(\"MIPscan_6pF\",\"low\",false,true,1,\"run_050575_injection_merged_ILCEnergyElectrons\"\)  

#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",true\) &
#root -l	-q analysis.cc\(\"3GeVMIPscan\",\"low\",true\) 

#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",false,true\) &
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",false,false\) &

#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,true\) 
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,false\)
