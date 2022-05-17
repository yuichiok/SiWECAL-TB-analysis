#sleep 20s
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",true\) &
#root -l	-q analysis.cc\(\"3GeVMIPscan\",\"low\",true\) 

root -l -q analysis.cc\(\"3GeVMIPscan_run_050060\",\"high\",false,true\)
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"high\",false,false\) &

#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,true\) 
#root -l -q analysis.cc\(\"3GeVMIPscan\",\"low\",false,false\)
