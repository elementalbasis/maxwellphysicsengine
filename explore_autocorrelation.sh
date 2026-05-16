#!/bin/zsh

for i in {000..009}; do
	julia demo-14.jl --particle-count 500 --particle-stiffness 10000 | tee "/home/WindowsShared/maxwell/run-2026-05-11_N=500_i=$i.tsv" > /dev/null
done
