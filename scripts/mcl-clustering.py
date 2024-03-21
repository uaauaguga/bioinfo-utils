#!/usr/bin/env python
import argparse
import math
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("MCL")
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(description='MCL clustering from pairwise search hits')
    parser.add_argument('--input', '-i', type=str, required=True, help='input filtered pairwise search hits')
    parser.add_argument('--output-directory','-od', type=str ,required=True, help="output directory")
    parser.add_argument('--metric','-m', type=str , choices=["bitscore","nle"], default="bitscore",help="similarity metric used for clustering")
    parser.add_argument('--threads', '-t', type=int, default=32, help='thread use for MCL clustering')
    parser.add_argument('--inflation', type=float, default=1.4, help='inflation parameter for MCL clustering')
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    
    logger.info("convert blast hits to abc format ...")    
    abc = os.path.join(args.output_directory,"pairwise.abc")
    fout = open(abc,"w")
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            query_id, hit_id = fields[:2]
            evalue, bitscore = float(fields[10]), int(fields[11])
            score = bitscore if args.metric == "bitscore" else max(-math.log10(evalue),200)
            print(query_id, hit_id, score, sep="\t", file=fout)
    fout.close()

    logger.info("load pairwise hits in abc format ...")
    mci =  os.path.join(args.output_directory,"pairwise.mci")
    tab = os.path.join(args.output_directory,"pairwise.tab")
    cmd = ["mcxload", "-abc", abc, "--stream-mirror", "-o", mci, "-write-tab", tab]
    subprocess.run(cmd)

    res = os.path.join(args.output_directory,"pairwise.mcl")
    logger.info("perform mcl clustering ...")
    cmd = ["mcl", mci, "-o", res, "-I", str(round(args.inflation,1)), "-te", str(args.threads)]
    subprocess.run(cmd)

    logger.info("dump clusters ...")
    clusters = os.path.join(args.output_directory,"pairwise.clstr")
    cmd = ["mcxdump", "-icl", res, "-tabr", tab, "-o", clusters]
    subprocess.run(cmd)

    logger.info("reforamt clusters ...")
    table = os.path.join(args.output_directory,"clusters.txt")
    fout = open(table,"w")
    with open(clusters) as fin:
        for i,line in enumerate(fin):
            fields = line.strip().split("\t")
            cluster_id = str(i).zfill(10)
            for seq_id in fields:
                print(seq_id, cluster_id, len(fields), sep="\t",file=fout)
    fout.close()
     
    
    
    

     
    

if __name__ == "__main__":
    main()
