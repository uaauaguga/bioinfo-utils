#!/usr/bin/env python
import argparse
from igraph import Graph
from collections import OrderedDict, defaultdict
import leidenalg as la
from tqdm import tqdm
import os
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('Leiden Community Detection')


def load_graph(path,max_evalue=0.01):
    node_index = 0
    node_index_lut = OrderedDict()
    edges = [] 
    with open(path) as f:
        for line in tqdm(f):
            # query, target, fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
            fields = line.strip().split("\t")
            query, target, fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue,bits = fields
            if float(evalue) > max_evalue:
                continue
            qstart, qend, tstart, tend = int(qstart), int(qend), int(tstart), int(tend)
            if query == target:
                continue
            if qstart > qend:
                #logger.warning(f"{query}-{target}: query start position is greater than end position, skip it .")
                continue
            if tstart > tend:
                #logger.warning(f"{query}-{target}: target start position is greater than end position, skip it .")
                continue
            if query in node_index_lut:
                index_1 = node_index_lut[query]
            else:
                index_1 = node_index
                node_index_lut[query] = index_1
                node_index += 1
            if target in node_index_lut:
                index_2 = node_index_lut[target]
            else:
                index_2 = node_index 
                node_index_lut[target] = index_2
                node_index += 1
            edges.append((index_1,index_2))
    g = Graph(edges) 
    g.vs["segments-ids"] = list(node_index_lut.keys())
    return g

def main():
    parser = argparse.ArgumentParser(description='Clustering pairwise search results with Leiden Algorithm')
    parser.add_argument('--input','-i',type=str,required=True,help="Input pairwise hits in blast m8 format")
    parser.add_argument('--output','-o',type=str,required=True, help="Where to save the membership")
    args = parser.parse_args()
    logger.info("Load pairwise hits into an igraph object ...")    
    g = load_graph(args.input)
    logger.info("Perform Leiden Partitioning ...")
    partition = la.find_partition(g, la.ModularityVertexPartition)
    g.vs["cluster"] = partition.membership
    logger.info("Saving membership assignment...")
    with open(args.output,"w") as f:
        f.write(f"segment_id\tdegree\tmembership\n")
        for segment_id, degree, membership in zip(g.vs["segments-ids"],g.vs.degree(),g.vs["cluster"]):
            f.write(f"{segment_id}\t{degree}\t{membership}\n")
    logger.info("All done .")

if __name__ == "__main__":
    main()
