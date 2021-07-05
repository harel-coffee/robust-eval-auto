

def run_robust(terminals, threshold, init, red):
    from amim_test_suite.algorithms.robust.pcst_approach.utils.ppi import PpiInstance, UnitEdgeWeight, read_ppi_shuffled
    from amim_test_suite.algorithms.robust.pcst_approach.utils import ExpMinMaxDiverseSteinerTreeComputer
    graph = read_ppi_shuffled("robustness_comparison/data/2020-07-07/protein-protein-interaction.txt", shuffle=True)
    edge_weights = UnitEdgeWeight()
    ppi_instance = PpiInstance(graph, terminals, edge_weights)
    engine = ExpMinMaxDiverseSteinerTreeComputer(initial_fraction=init,
                                                 reduction_factor=red)
    steiner_trees = engine(ppi_instance, n=10)
    t = steiner_trees.get_occurrences(include_terminals=True)
    if type(threshold) == list:
        return_dict = {}
        for thr in threshold:
            return_dict[thr] = t[t["%occurrences"] >= thr].index.values
        return return_dict
    return set(t[t["%occurrences"] >= threshold].index.values)