#!/usr/bin/env python
import os
import argparse
import sys
import subprocess
import numpy as np
import logging
from collections import defaultdict
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("motif discovery")

def main():
    parser = argparse.ArgumentParser(description='run motif discovery')
    parser.add_argument('--input',  '-i', type=str, required=True, help='input fasta')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='output location profile')
    parser.add_argument('--random-state', '-r', type=int, default=666,help='random seed')
    parser.add_argument('--min-width','-m',type=int,default=30,help="minimal width of the motif required")
    parser.add_argument('--max-width','-M',type=int,default=40,help="maximal width of the motif required")
    args = parser.parse_args()
    np.random.seed(args.random_state)
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    logger.info("load sequences ...")
    groupped_sequences = defaultdict(list)
    with open(args.input) as f:
        for header in f:
            sequence = next(f)
            if len(sequence.strip()) < args.max_width:
                continue
            rfam_id = header[1:].split("::")[0]
            groupped_sequences[rfam_id].append((header,sequence))
    rfam_ids = list(groupped_sequences.keys())
    n_train = int(len(rfam_ids)*0.8)
    train_ids = np.random.choice(rfam_ids,size=n_train,replace=False)
    n_validation = len(rfam_ids) - n_train
    validation_ids = [rfam_id for rfam_id in rfam_ids if rfam_id not in train_ids]
    logger.info(f"sample sequences for {n_train} families training ...")
    logger.info(f"the remained {n_validation} families are used for evaluation .") 
    train_fasta = os.path.join(args.output_directory, "train.fa")
    with open(train_fasta,"w") as f:
        for rfam_id in train_ids:
            i = np.random.randint(0,len(groupped_sequences[rfam_id]))
            header, sequence = groupped_sequences[rfam_id][i]
            f.write(header)
            f.write(sequence)
    validation_fasta = os.path.join(args.output_directory, "validation.fa")
    with open(validation_fasta,"w") as f:
        for rfam_id in validation_ids:
            i = np.random.randint(0,len(groupped_sequences[rfam_id]))
            header, sequence = groupped_sequences[rfam_id][i]
            f.write(header)
            f.write(sequence)
    motif = os.path.join(args.output_directory,"train","meme.txt") 
    if not os.path.exists(motif):
        logger.info("run motif finding ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/meme-env/bin/meme", "-mod", "zoops","-objfun", "classic", 
           "-oc", os.path.join(args.output_directory,"train"), "-rna", 
           "-minw", str(args.min_width), "-maxw", str(args.max_width), "-nmotifs", "5", train_fasta]
        subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    else:
        logger.info("results exist  ...")
    logger.info("convert motif ...")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/meme-env/bin/meme2meme", motif]
    minimal_motif = os.path.join(args.output_directory,"train","motif.meme") 
    fm = open(minimal_motif,"w")
    subprocess.run(cmd,stdout=fm)
    fm.close()   
    
    logger.info("filtering motifs ...")
    
    motifs = open(minimal_motif).read().split("\n\n\n")
    lines = motifs[0].split("\n")
    header = "\n".join(lines[:9])
    motifs[0] = "\n".join(lines[9:])
    minimal_motif_filtered = os.path.join(args.output_directory,"train","motif.filtered.meme")

    n = 0
    fm = open(minimal_motif_filtered, "w")
    fm.write(header+ "\n") 
    for motif in motifs:
        evalue = float(motif.split("\n")[9].split(" ")[-1]) 
        if evalue > 1:        
            continue
        n += 1
        fm.write(motif)
        fm.write("\n\n\n")
    fm.close()
    if n == 0:
        logger.info("No significant motif, exiting .")
        sys.exit(0)

    logger.info("select motif on validation set ...")
    logger.info("scan validation set ...")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/meme-env/bin/fimo","--norc",
           "--oc",os.path.join(args.output_directory, "validation"), minimal_motif_filtered , validation_fasta] 
    retcode = subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL) # stderr=subprocess.DEVNULL
    counter = defaultdict(int)
    hits = os.path.join(args.output_directory,"validation","fimo.tsv")

    seq_ids = set() 
    with open(hits) as f:
        _ = next(f) 
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            motif_id, seq_id = fields[1], fields[2]
            if seq_id in seq_ids:
                continue
            seq_ids.add(seq_id)
            counter[motif_id] += 1
    if len(counter) == 0:
        logger.info("no hit detected .")

    for motif_id in counter:
        print(motif_id, counter[motif_id], n_validation, counter[motif_id]/n_validation, sep="\t")
if __name__ == "__main__":
    main()
    
