#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=robustness2
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=250G

python robustness_comparison/robustness_mean_jaccard.py --algorithm MUST --threads 32
python robustness_comparison/robustness_mean_jaccard.py --algorithm RMUST --threshold 0.5 --trees 20 --threads 32

