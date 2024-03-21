#!/usr/bin/env python
import argparse
from collections import defaultdict
from tqdm import tqdm
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("groupping input sequence by a clustering table")

def main():
    parser = argparse.ArgumentParser(description='group sequence by clustering results')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequence in fasta format')
    parser.add_argument('--output','-o', type=str ,required=True, help="output groupped sequences")
    parser.add_argument('--representative','-r', type=str , help="output representative sequences")
    parser.add_argument('--table','-t', type=str ,required=True, help="clustering table")
    parser.add_argument('--min-size','-m', type=int , default = 0, help="minumn size of a cluster")
    parser.add_argument('--max-size','-M', type=int , default=10000000000, help="maximum size of a cluster")
    args = parser.parse_args()
    
    lookup = {}
    logger.info("load clustering table ...")
    with open(args.table) as f:
        for line in f:
            seq_id, rep_id = line.strip().split("\t")
            lookup[seq_id] = rep_id

    logger.info("groupping sequence ...")
    groupped_sequences = defaultdict(list)
    sequences = {}
    seq_id2header = {}
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                sequences[seq_id] = ""
                seq_id2header[seq_id] = line
            else:
                sequences[seq_id] += line.strip()
    for seq_id in sequences:
        if seq_id not in lookup:
            continue
        rep_id = lookup[seq_id]
        header = seq_id2header[seq_id]
        groupped_sequences[rep_id].append((header, sequences[seq_id] + "\n"))

    logger.info(f"{len(groupped_sequences)} clusters extracted .")
    logger.info("saving groupped sequences ...")
    n_too_small = 0
    n_too_large = 0
    n_passed = 0
    if args.representative is not None:
        frep = open(args.representative,"w")
    for i, rep_id in tqdm(enumerate(groupped_sequences)):
        index = str(i).zfill(8)
        fout = open(args.output + "/" + index + ".fa","w")
        if len(groupped_sequences[rep_id]) > args.max_size:
            n_too_large += 1
            continue
        if len(groupped_sequences[rep_id]) < args.min_size:
            n_too_small += 1
            continue
        n_passed += 1
        for header, sequence in groupped_sequences[rep_id]:
            if (args.representative is not None) and header[1:].startswith(rep_id):
                frep.write(header.strip() + " " + index + "\n")
                frep.write(sequence)
            fout.write(header.strip() + " " + index + "\n")
            fout.write(sequence)
        fout.close()
    if args.representative is not None:
        frep.close()
    logger.info(f"{n_passed} {n_too_small} {n_too_large}")
    logger.info("all done .")
        
    


if __name__ == "__main__":
    main()
