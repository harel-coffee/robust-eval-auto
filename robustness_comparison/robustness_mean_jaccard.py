import math
import os
from multiprocessing import Pool
import sys
import utils

from utils import get_parser, prep_parameters
from wrapper_robust import run_robust
from wrapper_must import run_must
from wrapper_rmust import run_rmust
from wrapper_diamond import run_diamond
from wrapper_domino import run_domino


def run_algorithm(algorithm, path_to_seeds, outfile, threshold=0.5, init=0.25, red=0.9):
    from amim_test_suite.algorithms.robust.pcst_approach.utils.ppi import read_terminals
    robustness_iterations = 20
    all_node_sets = [set() for i in range(robustness_iterations)]
    path_to_network = "robustness_comparison/data/2020-07-07/protein-protein-interaction.txt"
    if algorithm == 'ROBUST':
        terminals = read_terminals(path_to_seeds)
    for i in range(robustness_iterations):
        print(f'Iteration {i}')

        if algorithm == 'ROBUST':
            result_robust = run_robust(terminals, threshold, init, red)
            all_node_sets[i].update(result_robust)

        elif algorithm == 'DIAMOND':
            result_diamond = run_diamond(path_to_network, path_to_seeds, outfile)
            all_node_sets[i].update(result_diamond)

        elif algorithm == 'DOMINO':
            result_domino = run_domino(path_to_network, path_to_seeds, outfile)
            all_node_sets[i].update(result_domino)

        elif algorithm == 'MUST':
            result_must = run_must(path_to_network, path_to_seeds, outfile)
            all_node_sets[i].update(result_must)

        elif algorithm == 'RMUST':
            result_rmust = run_rmust(path_to_network, path_to_seeds, threshold)
            all_node_sets[i].update(result_rmust)

    #COMPUTE MEAN JACCARDS
    seed_name = path_to_seeds.split("/")[4].split(".")[0]
    sumJaccards = 0
    meanIntersectionSize = 0
    meanUnionSize = 0
    for i in range(robustness_iterations - 1):
        for j in range(i + 1, robustness_iterations):
            if len(all_node_sets[i].union(all_node_sets[j])) != 0:
                sumJaccards += (len(all_node_sets[i].intersection(all_node_sets[j])) / len(
                    all_node_sets[i].union(all_node_sets[j])))
                meanIntersectionSize += len(all_node_sets[i].intersection(all_node_sets[j]))
                meanUnionSize += len(all_node_sets[i].union(all_node_sets[j]))

    possibleCombinations = (1 / (math.factorial(robustness_iterations) // (
            math.factorial(2) * math.factorial(robustness_iterations - 2))))
    meanJaccard = possibleCombinations * sumJaccards
    meanIntersectionSize = possibleCombinations * meanIntersectionSize
    meanUnionSize = possibleCombinations * meanUnionSize
    if not os.path.exists(f'robustness_comparison/{algorithm}Out'):
        os.makedirs(f'robustness_comparison/{algorithm}Out')
    with open(outfile, "w") as file:
        file.write("Mean Jaccard: " + str(meanJaccard) + "\n")
        file.write("Mean Intersection Size: " + str(meanIntersectionSize) + "\n")
        file.write("Mean Union Size: " + str(meanUnionSize) + "\n")
        for nodeset in all_node_sets:
            file.write(str(nodeset) + "\n")
    if algorithm in ('ROBUST', 'RMUST'):
        return {f'{seed_name}_thr{threshold}': meanJaccard}
    else:
        return {seed_name: meanJaccard}


def call_algorithm(algorithm=utils.AlgorithmSelector.ROBUST, threshold=0.5, init=0.25, red=0.9, threads=4):
    print(f'Running {algorithm} in {threads} threads...')
    list_tuples = prep_parameters(algorithm, threshold, init, red)
    all_jaccards = {}

    with Pool(threads) as p:
        return_dict = p.starmap(run_algorithm, list_tuples)

    for d in return_dict:
        all_jaccards.update(d)

    with open(f'robustness_comparison/{algorithm.value}Out/{algorithm}.out', "w") as file:
        file.write("seed_set,mean jaccard\n")
        for key, value in all_jaccards.items():
            file.write(f'{key},{value}\n')
        file.close()


if __name__ == '__main__':
    sys.path.append(os.getcwd())
    arguments = get_parser().parse_args()
    if not arguments.threshold:
        threshold = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    else:
        threshold = arguments.threshold
    call_algorithm(arguments.algorithm, threshold, arguments.init, arguments.red, arguments.threads)