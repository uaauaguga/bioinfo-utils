#!/usr/bin/env python
from ete3 import Tree
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('rename tree node')

def main():
    parser = argparse.ArgumentParser(description='rename tree node')
    parser.add_argument('--input','-i', type=str, required=True, help='input tree in newick format')
    parser.add_argument('--name-lookup','-nl',type=str,required=True, help="name lookup")
    parser.add_argument('--output','-o',type=str, required=True, help = 'output tree')
    parser.add_argument('--format','-f',type=int, default=1, help = 'format defined by ete package')
    args = parser.parse_args() 

    tree = Tree(newick=args.input,format=args.format,quoted_node_names=True)
    rename = {}
    with open(args.name_lookup) as f:
        for line in f:
            name1, name2 = line.strip().split("\t")
            rename[name1] = name2
    for node in tree.traverse():
        if node.name in rename:
            node.name = rename[node.name]
    with open(args.output,"w") as f:
        f.write(tree.write(format=args.format))



if __name__ == "__main__":
    main()
