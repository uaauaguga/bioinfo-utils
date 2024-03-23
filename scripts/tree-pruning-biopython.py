#!/usr/bin/env python
from Bio import Phylo
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('tree pruning')
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='pruning a tree in newick format')
    parser.add_argument('--input','-i',type=str, required=True, help='input tree in newick format')
    parser.add_argument('--nodes','-n',type=str, required = True, help="name of nodes to extract")
    parser.add_argument('--output','-o',type=str, required=True, help = 'output tree')
    args = parser.parse_args() 

    logger.info("load tree ...")
    tree = Phylo.read(args.input, "newick")

    logger.info("load node names ...")
    node_ids = set(open(args.nodes).read().strip().split("\n"))
    nodes_to_prune = []

    logger.info("extract node to prune ...")
    for clade in tree.get_terminals():
        if clade.name not in node_ids:
            nodes_to_prune.append(clade)

    logger.info("pruning the tree ...")
    for node in tqdm(nodes_to_prune):
        tree.prune(node)

    logger.info("Saving results ...")
    with open(args.output,"w") as fout:
        fout.write(tree.format("newick"))

    logger.info("All done .")

if __name__ == "__main__":
    main()
