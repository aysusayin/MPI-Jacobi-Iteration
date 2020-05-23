#!/usr/bin/env bash

mpicxx main.cpp -o mpi_project

echo "Program compiled"

lines=('N,T1,T8,T27,T64,S8,S27,S64')
tmp_file_name="tmp_file"

for N in 11 35 59 71; do
    for t in 1 8 27 64; do
          /usr/bin/time -o $tmp_file_name -f "%e" mpiexec -n $t ./mpi_project --silent $N
          exec_times[$t]="$(grep -oP "\d+.*\d+" $tmp_file_name)"
          echo "Analysis for N=$N with $t processes is done."
          
    done
    s8=$(bc <<<"scale=3; ${exec_times[1]} / ${exec_times[8]}")
    s27=$(bc <<<"scale=3; ${exec_times[1]} / ${exec_times[27]}")
    s64=$(bc <<<"scale=3; ${exec_times[1]} / ${exec_times[64]}")
    line="$N,${exec_times[1]},${exec_times[8]},${exec_times[27]},${exec_times[64]},$s8,$s27,$s64"
    lines+=("${line}")
  
done

file_name="MPI_Analysis.csv"

rm -f ${tmp_file_name}
rm -f ${file_name}

for line in "${lines[@]}"; do
  echo "$line" >>"$file_name"
done

echo "Done!"
