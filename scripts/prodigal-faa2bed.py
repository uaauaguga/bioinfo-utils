#!/usr/bin/env python
import argparse
import gzip
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract bed from aa file')

def main():
    parser = argparse.ArgumentParser(description='protein sequence to bed')
    parser.add_argument('--input', '-i', required=True, help="input protein sequences")
    parser.add_argument('--output','-o', required=True, help="output intervals in bed format")
    args = parser.parse_args() 

    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            # >NZ_JXRQ01000002.1_2 # 1377 # 2840 # 1 # ID=2_2;partial=00;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10bp;gc_cont=0.460
            fields = line[1:].strip().split(" ")
            protein_id, start, end, strand = fields[0], fields[2], fields[4], fields[6]
            strand = {"1":"+","-1":"-"}[strand]
            contig_id = protein_id[:protein_id.rfind("_")]
            start, end = int(start)-1, int(end)
            print(contig_id, start, end, protein_id, ".", strand,file=fout,sep="\t")
    fout.close()
           

if __name__ == "__main__":
    main()
