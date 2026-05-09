#!/bin/zsh

a=$1
b=$2
N=500
for k in 10000 15000 20000 25000 30000 35000 40000 45000 50000 55000 60000; do
	for i in {$a..$b}; do
		julia demo-13.jl --particle-count $N --particle-stiffness=$k | tee "/home/WindowsShared/maxwell/run-2026-05-09_N=${N}_k=${k}_i=${i}.tsv" > /dev/null
	done
done
