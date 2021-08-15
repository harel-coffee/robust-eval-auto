#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --job-name=amim
#SBATCH --output=%x.%j.txt
#SBATCH --error=%x.%j.err
#SBATCH --mem=200G

python amim_test_suite/run_tests.py parallel --networks BIOGRID HPRD STRING APID IID_BRAIN IID_LUNG --generators ORIGINAL REWIRED EXPECTED_DEGREE SHUFFLED SCALE_FREE UNIFORM --method MUST --verbose
