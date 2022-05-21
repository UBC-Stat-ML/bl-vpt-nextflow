#!/bin/bash

module load git
module load singularity
module load miniconda3

./nextflow -c pbspro.config $@ -with-report -with-trace -with-timeline -with-dag -qs 1000
