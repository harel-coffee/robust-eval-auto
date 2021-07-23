#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --job-name=diamondrobustness
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=200G

python robustness_comparison/robustness_mean_jaccard.py --algorithm DIAMOND --threads 25

