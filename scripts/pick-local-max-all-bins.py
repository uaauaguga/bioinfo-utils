#!/usr/bin/env python
import argparse
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("select intervals")

from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='pick best interval from overlapped ones')
    parser.add_argument('--input', '-i',required=True,help="input intervals")
    parser.add_argument('--output','-o',required=True,help="output intervals")
    args = parser.parse_args()


    logger.info(f"load intervals from {args.input} ...")
    logger.info(f"picked intervals will be saved to {args.output} .")
    fin = open(args.input)
    fout = open(args.output,"w")    
    starts, ends = defaultdict(list), defaultdict(list)
    strands = defaultdict(list)
    scores = defaultdict(list)
    names = defaultdict(list)
    for line in fin:
        fields = line.strip().split("\t")
        seq_id, start, end, name, score, strand = fields[:6]
        start, end, score = int(start), int(end), float(score)
        starts[seq_id].append(start)
        ends[seq_id].append(end)
        strands[seq_id].append(strand)
        scores[seq_id].append(score)
        names[seq_id].append(name)
    #lengths = {} 
    #for seq_id in ends:
    #    lengths[seq_id] = max(ends[seq_id]) 
    
    for seq_id in scores:
        values = np.array(scores[seq_id])
        shift_right = np.array([0] + scores[seq_id])
        shift_left = np.array(scores[seq_id] + [0])
        mask = (values > shift_right[:-1]) & (values > shift_left[1:])
        print(mask.sum())
        indices = np.where(mask)[0]
        for i in indices:
            start, end, name, score, strand = starts[seq_id][i], ends[seq_id][i], names[seq_id][i], scores[seq_id][i], strands[seq_id][i]
            print(seq_id, start, end,name,score,strand,sep="\t",file=fout)
    fin.close()
    fout.close()

if __name__ == "__main__":
    main()

