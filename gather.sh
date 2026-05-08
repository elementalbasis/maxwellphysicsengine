#!/bin/bash

for k in 50000 30000 10000 5000 3000 1000 500 300 100; do
	julia demo-12.jl --particle-stiffness=$k > "/home/WindowsShared/maxwell/run-2026-05-07_k=$k.tsv"
done
