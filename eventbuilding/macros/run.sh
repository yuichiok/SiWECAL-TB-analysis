run=32001
root -l -q number_occurencies.C\(\"$run\"\)
root -l -q time_correlation.C\(\"$run\"\)

run=32002
root -l -q number_occurencies.C\(\"$run\"\)
root -l -q time_correlation.C\(\"$run\"\)

run=32004
root -l -q number_occurencies.C\(\"$run\"\)
root -l -q time_correlation.C\(\"$run\"\)

run=32015
root -l -q number_occurencies.C\(\"$run\"\)
root -l -q time_correlation.C\(\"$run\"\)

root -l -q time_correlation_all_files.C
