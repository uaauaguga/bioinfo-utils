#!/usr/bin/env python
from ete3 import Tree
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('tree pruning')

def main():
    parser = argparse.ArgumentParser(description='pruning a tree in newick format')
    parser.add_argument('--input','-i',type=str, required=True, help='input tree in newick format')
    parser.add_argument('--nodes','-n',type=str, required = True, help="name of nodes to extract")
    parser.add_argument('--output','-o',type=str, required=True, help = 'output tree')
    parser.add_argument('--format','-f',type=int, default=1, help = 'format defined by ete package')
    args = parser.parse_args() 
    logger.info("load tree ...")
    tree = Tree(newick=args.input,format=args.format,quoted_node_names=True)
    logger.info("load node names ...")
    names = set(open(args.nodes).read().strip().split("\n"))
    nodes = []
    logger.info("look for specified nodes ...")
    for node in tree.traverse():
        name = node.name
        if name in names:
            nodes.append(node)
    logger.info(f"{len(nodes) in {len(names)}} specified nodes present in the tree .")
    logger.info("pruning the tree ...")
    tree.prune(nodes)
    logger.info("saving output ...")
    with open(args.output,"w") as f:
        f.write(tree.write(format=args.format))
    logger.info("all done .")



if __name__ == "__main__":
    main()
