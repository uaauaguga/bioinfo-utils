#!/usr/bin/env python
import argparse
import re
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Reformat cd-hit clustering table')
    parser.add_argument('--input', '-i', required=True, help="Input path")
    parser.add_argument('--output','-o', required=True, help="Output path")
    parser.add_argument('--type','-t', required=True, choices = ["aa","nt"],help="Data type, either nucleotide or amino acid")
    args = parser.parse_args()
    fin = open(args.input)
    fout = open(args.output,"w")
    if args.type == "nt":
        pattern = r"(\d+)nt, >(.+)\.\.\. (.+)" 
    else:
        pattern = r"(\d+)aa, >(.+)\.\.\. (.+)" 
    seq_ids = []
    centroid_id = ""
    for line in tqdm(fin):
        if line.startswith(">"):
            if len(seq_ids) > 0:
                assert len(centroid_id) > 0
                for seq_id in seq_ids:
                    fout.write(f"{seq_id}\t{centroid_id}\n")
            seq_ids = []
            centroid_id = "" 
        else:
            cluster_id, description  = line.strip().split("\t")
            m = re.match(pattern, description) 
            length, seq_id, centroid = m.groups()
            seq_ids.append(seq_id)
            if centroid == "*":
                centroid_id = seq_id
    if len(seq_ids) > 0:
        assert len(centroid_id) > 0
        for seq_id in seq_ids:
            fout.write(f"{seq_id}\t{centroid_id}\n")
    fin.close()
    fout.close()
if __name__ == "__main__":
    main()    
