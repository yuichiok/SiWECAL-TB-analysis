#sleep 20s
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",true\) &
#root -l	-q analysis.cc\(\"3GeVMIPscan\",\"low\",true\) 

# root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",false,true\) &
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",false,false\) &

#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,true\) 
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,false\)

root -l -q analysis.cc\(\"Run_ILC_03112022_cosmic_it16\",\"high\",false,true\)