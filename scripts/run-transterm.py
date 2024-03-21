#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
import re
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("run transterm-HP")
import io


def main():
    parser = argparse.ArgumentParser(description='run transterm-HP for terminator prediction')
    parser.add_argument('--input','-i',type=str, required=True, help='input fasta file')
    parser.add_argument('--output','-o',type=str, required=True, help="output directory")
    parser.add_argument('--selection', '-s', type=str, default="local",
                        choices=["local","full"], help = "select locla best scoring interval or best interval of the whole sequence")
    args = parser.parse_args()

    logger.info(f"predict terminator for {args.input} ...")
    logger.info(f"results will be saved to {args.output} .")

    if not os.path.exists(args.output):
        logger.info(f"{args.output} does not exists, create it .")
        os.makedirs(args.output)

    logger.info("make a fake coordinate for transterm-HP ...")
    fin = open(args.input)
    sequences = {} 
    for line in fin:
        if line.startswith(">"):
            attrs = line.strip()[1:].split(" ")
            seq_id = attrs[0]
            #assert seq_id not in sequences, "sequence with identical id detected ."
            sequences[seq_id] = ""
        else:
            sequences[seq_id] += line.strip()
            
    fcrd = open(os.path.join(args.output,"fake.crd"),"w")
    fseq = open(os.path.join(args.output,"sequence.gt74.fa"),"w")
    n_too_short = 0
    for seq_id in sequences:
        sequence = sequences[seq_id]
        L = len(sequence)
        if L < 75:
            n_too_short += 1
            continue
        fseq.write(f">{seq_id}\n")
        fseq.write(f"{sequence}\n")
        name = f"{seq_id}-fake-1"
        fcrd.write(f"{name}\t1\t2\t{seq_id}\n")
        name = f"{seq_id}-fake-2"
        fcrd.write(f"{name}\t{3}\t{L-2}\t{seq_id}\n")
        name = f"{seq_id}-fake-3"
        fcrd.write(f"{name}\t{L-1}\t{L}\t{seq_id}\n")
    fseq.close()
    fcrd.close()

    logger.info(f"{n_too_short} sequences shorter than 75 nt are discarded .")

    logger.info("run transterm-HP ...")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/tools/transterm_hp_v2.09/transterm","--min-conf=50","--all-context",
           "-p","/apps/home/lulab_jinyunfan/qhsky1/tools/transterm_hp_v2.09/expterm.dat",
           os.path.join(args.output,"sequence.gt74.fa"),
           os.path.join(args.output,"fake.crd")]
    fout = open(os.path.join(args.output,"transterm.txt"),"w")
    logger.info(" ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
    for line in io.TextIOWrapper(proc.stdout, encoding="unicode_escape"):
        fout.write(line)
    fout.close() 
    code =  proc.poll()
    assert code == 0, f"transterm failed with {code} ."
    logger.info("convert results to bed format ...")
    bed = os.path.join(args.output,"transterm.bed")
    fout = open(bed,"w")
    with open(os.path.join(args.output,"transterm.txt"),encoding = "ISO-8859-1") as f:
        for line in f:
            if line.startswith("SEQUENCE"):
                line = line.strip()
                contig_id = re.split("\s+",line)[1]
            elif line.startswith("  TERM"):
                line = line.strip() 
                fields = re.split("\s+",line)
                start, end = int(fields[2]), int(fields[4])
                start -= 1
                strand = fields[5]
                if start > end:
                    end, start = start, end
                    assert strand == "-"
                score = fields[7]
                print(contig_id,start,end,".",score,strand,sep="\t",file=fout)
    fout.close()
    logger.info("sort bed file ...")   
 
    proc = subprocess.run(["sort","-k1,1","-k2,2n","-o",bed,bed])

    if args.selection == "local":
        logger.info("merge overlapped intervals, and assign the best score ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/python","scripts/pick-local-max.py", "-i", os.path.join(args.output,"transterm.bed"),
               "-o",os.path.join(args.output,"transterm.max.bed")]
        logger.info(" ".join(cmd))
        proc = subprocess.run(cmd)
        #code =  proc.poll()
        #assert code == 0, f"merge intervals failed with {code} ."
    else:
        logger.info("take highes scoring interval from each sequence ...")
        max_score_by_sequence = defaultdict(int)
        ivs = {}
        with open(os.path.join(args.output,"transterm.bed")) as f:
            for line in f:
                fields = line.strip().split("\t")
                seq_id = fields[0]
                score = int(fields[4])
                if score > max_score_by_sequence[seq_id]:
                    ivs[seq_id] = fields[1:]
                    max_score_by_sequence[seq_id] = score
        fout = open(os.path.join(args.output,"transterm.max.bed"),"w")    
        for seq_id in ivs:
            print(seq_id,*ivs[seq_id],sep="\t",file=fout)
        fout.close()
    logger.info("All done .")


if __name__ == "__main__":
    main() 
