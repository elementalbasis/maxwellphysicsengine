#!/bin/zsh

a=$1
b=$2
N=500
#for k in 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000; do
#for k in 65000 70000 75000 80000 85000 90000 100000; do
#for k in 3000000 5000000 10000000; do
#	for i in {$a..$b}; do
#		#julia demo-13.jl --particle-count $N --particle-stiffness=$k | tee "/home/WindowsShared/maxwell/run-2026-05-09_N=${N}_k=${k}_i=${i}.tsv" > /dev/null
#		julia demo-13.jl --particle-count $N --particle-stiffness=$k --simulation-dt "5e-7" | tee "/home/WindowsShared/maxwell/run-2026-05-09_N=${N}_k=${k}_i=${i}.tsv" > /dev/null
#	done
#done
#
for m in 0.00005 0.0001 0.0003 0.0005 0.001 0.003 0.005 0.01 0.03 0.05 0.1 0.3 0.5; do
	for i in {$a..$b}; do
		julia demo-13.jl --particle-count $N --particle-mass $m > "/home/WindowsShared/maxwell/run-2026-05-15_N=${N}_m=${m}_i=${i}.tsv"
	done
done
