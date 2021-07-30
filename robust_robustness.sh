#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=robust_robustness
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=100G

for initial_fraction in 0.25 0.5 0.75
do
  for reduction_factor in 0.1 0.3 0.5 0.7 0.9
  do
    python robustness_comparison/robustness_mean_jaccard.py --algorithm ROBUST --threshold 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 --init $initial_fraction --red $reduction_factor --threads 10 --trees 5 10 15 20
  done
done
