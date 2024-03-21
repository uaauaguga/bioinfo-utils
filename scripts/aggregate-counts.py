#!/usr/bin/env python
from collections import defaultdict
import argparse

def main():
    parser = argparse.ArgumentParser(description='aggregate counts by a specified sequence groupping')
    parser.add_argument('--input', '-i', type=str, required = True,help="input counts")
    parser.add_argument('--output','-o', type=str, required = True, help="groupped counts")
    parser.add_argument('--table', '-t', type=str, 
      default = "/apps/home/lulab_jinyunfan/qhsky1/sRNA-finding/output/homolog-search/combined-stranded/candidate-tx-level3-clustering-table.txt", help="groupping table")
    parser.add_argument("--black-list","-bl", type = str,
      default="/apps/home/lulab_jinyunfan/qhsky1/sRNA-finding/output/homolog-search/combined/candidate-tx-groupped-by-level3.nodupe.lowcomplexcity.fa",help="black list to exclude")
    args = parser.parse_args()

    counts_by_sequence_fwd = {}
    counts_by_sequence_rev = {}
    with open(args.input) as f:
        for line in f:
            seq_id, count_fwd, count_rev = line.strip().split("\t")
            counts_by_sequence_fwd[seq_id] = int(count_fwd)
            counts_by_sequence_rev[seq_id] = int(count_rev) 

    black_list = set()
    with open(args.black_list) as f:
        for line in f:
            seq_id = line.strip().split("\t")[0]
            black_list.add(seq_id)

    counts_by_cluster_fwd = defaultdict(int)    
    counts_by_cluster_rev = defaultdict(int)
    cluster_ids = set()
    with open(args.table) as f:
        for line in f:
            seq_id, cluster_id = line.strip().split("\t")
            if seq_id not in counts_by_sequence_fwd:
                continue
            if seq_id in black_list:
                continue
            cluster_ids.add(cluster_id)
            counts_by_cluster_fwd[cluster_id] += counts_by_sequence_fwd[seq_id]    
            counts_by_cluster_rev[cluster_id] += counts_by_sequence_rev[seq_id]
    
    fout = open(args.output,"w")
    for cluster_id in sorted(list(cluster_ids)):
        print(cluster_id, counts_by_cluster_fwd[cluster_id], counts_by_cluster_rev[cluster_id], sep="\t",file=fout)
    fout.close()
        
    


if __name__ == "__main__":
    main()
