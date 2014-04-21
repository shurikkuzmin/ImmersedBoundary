#!/bin/bash
rm vtk*.vtk
rm vel*.dat
g++ -O3 grid.cpp -o main.out
