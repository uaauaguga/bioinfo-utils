#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
import logging
import io
from utils import merge, subtract
import pyfaidx
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('reduce redundancy')

def main():
    parser = argparse.ArgumentParser(description='remove nearly identical sequence by minimap2 search')
    parser.add_argument('--input', '-i', required=True, help="input sequences")
    parser.add_argument('--output','-o', required=True, help="output sequences")
    parser.add_argument('--min-length','-ml', type=int, default=1000, help="min length of a genome segment to consider")
    parser.add_argument('--thread','-t',type=int,default=3,help="thread for minimap2 mapping")
    args = parser.parse_args()

    if not os.path.exists(args.input + ".fai"):
        logger.info("index input sequence ...")
        cmd = ["samtools","faidx",args.input]
        subprocess.run(cmd)

    lengths = {}
    with open(args.input + ".fai") as f:
        for line in f:
            seq_id ,length = line.strip().split("\t")[:2]
            length = int(length)
            lengths[seq_id] = length

    logger.info("run minimap2 for pairwise search")
    cmd = ["minimap2","-X","-t",str(args.thread),"-x","asm20",args.input, args.input]
    print(" ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
    handle = io.TextIOWrapper(proc.stdout, encoding="utf-8")
    
    intervals = {}
    for line in handle:
        fields = line.strip().split("\t")
        qname, tname = fields[0], fields[5]
        qlength, tlength = int(fields[1]), int(fields[6])
        # query interval field 3, 4
        qstart, qend = int(fields[2]), int(fields[3])
        # target interval field 8, 9
        tstart, tend = int(fields[7]), int(fields[8])
        aligned_length = int(fields[9]) # define aligned length as number of residue matches
        divergence = None
        for field in fields[12:]:
            if field.startswith("dv"):
                divergence = float(field.split(":")[-1])
                break
        assert divergence is not None, f"{line.strip()} does not contain a dv tag ..."            
        if divergence > 0.01:
            continue 
        # always remove sequence with larger id
        if qname > tname:
            seq_id = qname
            start, end = qstart, qend
        else:
            seq_id = tname
            start, end = tstart, tend
        if seq_id not in intervals:
            intervals[seq_id] = []
        # apend to intervals to considered as redundant
        intervals[seq_id].append((start, end))
    full_length = 0
    kept_length = 0 
    fasta = pyfaidx.Fasta(args.input)
    logger.info("extract sequences ...")
    fout = open(args.output,"w")
    for seq_id in lengths:
        full_length += lengths[seq_id]
        if seq_id not in intervals:
            if lengths[seq_id] < args.min_length:
                continue
            sequence = str(fasta[seq_id][::])
            fout.write(f">{seq_id}:0-{lengths[seq_id]}\n")
            p = 0
            while p < len(sequence):
                fout.write(sequence[p:p+70] + "\n")
                p += 70
            kept_length += lengths[seq_id]
        else:
            # exlude sequence considered as redundant
            ivs = merge(intervals[seq_id])
            ivs = subtract(ivs, lengths[seq_id])      
            for start, end in ivs:
                if end - start < args.min_length:
                    continue
                kept_length += end - start
                sequence = str(fasta[seq_id][start:end])
                fout.write(f">{seq_id}:{start}-{end}\n")
                p = 0
                while p < len(sequence):
                    fout.write(sequence[p:p+70] + "\n")
                    p += 70
    logger.info(f"{full_length} {kept_length} {round(kept_length/full_length,4)}")
    fout.close()
    logger.info("All done.")

if __name__ == "__main__":
    main()
