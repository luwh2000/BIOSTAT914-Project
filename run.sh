#!/bin/bash

input=${1:-model}

file_number=${2:-1}
./$input sample num_samples=1000 num_warmup=1000 data file=input/gene_${file_number}.json output file=output/output_${file_number}.csv
