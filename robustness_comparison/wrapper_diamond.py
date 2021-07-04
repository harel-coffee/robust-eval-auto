import os
import pandas as pd
import random
from subprocess import call


def run_diamond(path_to_network, path_to_seeds, outfile):
    seed_name = path_to_seeds.split("/")[4].split(".")[0]
    file = open(path_to_network)
    lines = file.readlines()
    lines.pop(0)
    random.shuffle(lines)
    tmpfile = f'shuffled_network{seed_name}.txt'
    f = open(tmpfile, "w")
    f.writelines(lines)
    f.close()
    if not os.path.exists('DiamondOut'):
        os.makedirs('DiamondOut')
    call(["python3", "DIAMOnD.py", tmpfile, path_to_seeds, "200", "1", outfile])
    nodes = pd.read_csv(outfile, sep="\t")["DIAMOnD_node"]
    os.remove(tmpfile)
    os.remove(outfile)
    return nodes