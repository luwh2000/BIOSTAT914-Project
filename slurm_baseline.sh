#!/bin/bash

for i in {1..1000}
do
    sbatch --wrap="source ./run.sh model_baseline $i"
done
