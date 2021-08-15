#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=runtime
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=200G

python runtime_evaluation/compare_runtimes.py