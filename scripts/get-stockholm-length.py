#!/usr/bin/env python
import argparse
import os
from Bio import AlignIO
import numpy as np
from tqdm import tqdm

def get_length(path):
    align = AlignIO.read(path, "stockholm")
    lengths = []
    for record in align:
        sequence = str(record.seq).replace("-","").replace(".","")
        lengths.append(len(sequence))
    mean = int(np.mean(lengths))
    median = int(np.median(lengths))
    return mean, median, len(lengths)

def main():
    parser = argparse.ArgumentParser(description='Get length of alignments in stockholm format')
    parser.add_argument('--indir', '-i',required=True,help="Input dir contains stk files")
    parser.add_argument('--output','-o',required=True,help="Output table of average length")
    args = parser.parse_args()
    fout = open(args.output,"w")
    fout.write("rfam-id\tmean\tmedian\tnumber\n")
    for stk in tqdm(os.listdir(args.indir)):
        rfam_id = stk.split(".")[0]
        path = os.path.join(args.indir,stk)
        mean, median, number  = get_length(path)
        print(rfam_id, mean, median, number, sep="\t",file=fout)
    fout.close()
    
   
    
if __name__ == "__main__":
    main()

