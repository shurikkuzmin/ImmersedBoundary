#!/bin/bash
rm vtk*.vtk
rm vel*.dat
#g++ -O3 moving_cylinder.cpp -o main.out
g++ -O3 ellipse.cpp -o main.out
