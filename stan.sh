#!/bin/bash

# Set the first input parameter with a default value of mode.stan
input=${1:-model}

# Activate the conda environment
cd /hpc/home/wl287/BIOSTAT914/biostat914/bin/cmdstan
make /hpc/home/wl287/BIOSTAT914/$input
cd /hpc/home/wl287/BIOSTAT914
