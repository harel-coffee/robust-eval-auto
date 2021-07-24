import os
from enum import Enum
import argparse


class AlgorithmSelector(Enum):
    ROBUST = 'ROBUST'
    DOMINO = 'DOMINO'
    DIAMOND = 'DIAMOND'
    MUST = 'MUST'
    RMUST = 'RMUST'

    def __str__(self):
        return self.value


def allowed_values_floats(x):
    x = float(x)
    if not 0.0 <= x <= 1.0:
        raise argparse.ArgumentTypeError("Choose a value between 0 and 1")
    return x


def get_parser():
    parser = argparse.ArgumentParser('tests the mean overlap between different runs with shuffled input networks')
    parser.add_argument('--algorithm', type=AlgorithmSelector, choices=list(AlgorithmSelector), required=True)
    parser.add_argument('--threshold', type=allowed_values_floats, nargs='+',
                        help='Threshold for ROBUST and R-MuST, default: 0.1-0.9 in 0.1 steps')
    parser.add_argument('--init', default=0.25, type=allowed_values_floats,
                        help='Initial fraction value for ROBUST, default: 0.25')
    parser.add_argument('--red', default=0.9, type=allowed_values_floats,
                        help='Reduction factor value for ROBUST, default: 0.9')
    parser.add_argument('--threads', default=4, type=int, help='Run on multiple seed files at once with x threads')
    parser.add_argument('--trees', default=10, type=int, help='Number of trees to be returned.', nargs='+')
    return parser


def prep_parameters(algorithm=AlgorithmSelector.ROBUST, threshold=0.5, init=0.25, red=0.9, nr_of_trees=20):
    list_tuples = []
    #sort seeds by size
    path = "robustness_comparison/data/2020-07-07/all-seeds"
    dirpath = os.path.abspath(path)
    all_files = (os.path.join(basedir, filename) for basedir, dirs, files in os.walk(dirpath) for filename in files)
    all_seeds = sorted(all_files, key = os.path.getsize)
    #tmp
    all_seeds = [os.path.basename(file) for file in all_seeds][1:900]
    if type(threshold) == list and len(threshold) == 1:
        threshold = threshold[0]
    if type(nr_of_trees) == list and len(nr_of_trees) == 1:
        nr_of_trees = nr_of_trees[0]

    if type(threshold) == list and algorithm in (AlgorithmSelector.ROBUST, AlgorithmSelector.RMUST):
        print(f'Testing {algorithm} with multiple thresholds...')
        if type(nr_of_trees) == list:
            print(f'Testing {algorithm} with multiple numbers of trees...')
            for seed_file in all_seeds:
                path_to_seeds = f"robustness_comparison/data/2020-07-07/all-seeds/{seed_file}"
                if algorithm == AlgorithmSelector.ROBUST:
                    robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_init{init}_red{red}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red, nr_of_trees))
                else:
                    raise ValueError("#trees can only be specified as a list for ROBUST")
        else:
            for seed_file in all_seeds:
                path_to_seeds = f"robustness_comparison/data/2020-07-07/all-seeds/{seed_file}"
                if algorithm == AlgorithmSelector.ROBUST:
                    robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_init{init}_red{red}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red))
                else:
                    rmust_out = f"robustness_comparison/{algorithm}Out/RMUST_{seed_file.split('.')[0]}.out"
                    list_tuples.append((algorithm.value, path_to_seeds, rmust_out, threshold, nr_of_trees))
    elif algorithm == AlgorithmSelector.ROBUST:
        print(f'Running ROBUST with initial fraction={init} and reduction factor={red}')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/2020-07-07/all-seeds/{seed_file}"
            robust_out = f"robustness_comparison/{algorithm}Out/ROBUST_{seed_file.split('.')[0]}_thr{threshold}_init{init}_red{red}.out"
            list_tuples.append((algorithm.value, path_to_seeds, robust_out, threshold, init, red, nr_of_trees))
    elif algorithm in [AlgorithmSelector.DIAMOND, AlgorithmSelector.DOMINO, AlgorithmSelector.MUST]:
        print(f'Running {algorithm}...')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/2020-07-07/all-seeds/{seed_file}"
            outfile = f"robustness_comparison/{algorithm}Out/{algorithm}_{seed_file.split('.')[0]}.out"
            list_tuples.append((algorithm.value, path_to_seeds, outfile))
    elif algorithm == AlgorithmSelector.RMUST:
        print(f'Running R-MuST with threshold={threshold}')
        for seed_file in all_seeds:
            path_to_seeds = f"robustness_comparison/data/2020-07-07/all-seeds/{seed_file}"
            rmust_out = f"robustness_comparison/{algorithm}Out/RMUST_{seed_file.split('.')[0]}_thr{threshold}.out"
            list_tuples.append((algorithm.value, path_to_seeds, rmust_out, threshold, nr_of_trees))
    return list_tuples
