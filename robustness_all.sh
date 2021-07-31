#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --job-name=robustness1
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=250G

python robustness_comparison/robustness_mean_jaccard.py --algorithm ROBUST --threshold 0.1 --init 0.25 --red 0.9 --threads 32 --trees 30
python robustness_comparison/robustness_mean_jaccard.py --algorithm RMUST --threshold 0.1 --trees 30 --threads 32
#python robustness_comparison/robustness_mean_jaccard.py --algorithm MUST --threads 32
#python robustness_comparison/robustness_mean_jaccard.py --algorithm DIAMOND --threads 25
#python robustness_comparison/robustness_mean_jaccard.py --algorithm DOMINO --threads 25

