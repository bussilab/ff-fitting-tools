#!/bin/bash

lepton_path="./include/lepton"

make clean

#Compile Lepton
(
cd $lepton_path ;
	g++ -c *.cpp -I.. -fopenmp
)

make
