#!/bin/bash

module load gcc/5.4.0
module load git
module load singularity

./nextflow -c sockeye.config $@