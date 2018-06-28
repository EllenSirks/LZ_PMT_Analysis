#!/bin/bash

cd /home/ellen/Desktop/MPhys_Project/Output/

if [ "$1" = "" ]
then 
	echo No input folder given.
	exit
fi

echo You entered the folder: $1

mkdir $1'_output'
cd $1'_output'
mkdir A
mkdir B

cd A
mkdir histograms
mkdir specs
mkdir plots

cd plots
mkdir count
mkdir chi2

cd ../../B
mkdir histograms
mkdir specs
mkdir plots

cd plots
mkdir count
mkdir chi2