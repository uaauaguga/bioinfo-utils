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
logger = logging.getLogger('extract sequences')

def main():
    parser = argparse.ArgumentParser(description='extract sequence not present in reprepsentative genome')
    parser.add_argument('--input-directory', '-i', required=True, help="input directory containing genomes of s species")
    parser.add_argument('--repsentative-genome', '-r', required=True, help="representative genome, should present in input directory")
    parser.add_argument('--output','-o', required=True, help="output path")
    parser.add_argument('--min-length','-ml', type=int, default=1000, help="min length of a genome segment to consider")
    args = parser.parse_args()

    rep_genome_path = os.path.join(args.input_directory, args.repsentative_genome)
    if not os.path.exists(rep_genome_path):
        logger.error(f"{rep_genome_path} does not exists, exiting .")
    fout = open(args.output,"w")
    for genome in os.listdir(args.input_directory):
        if not genome.endswith(".fna"):
            continue
        intervals = {}
        if genome == args.repsentative_genome:
            continue
        logger.info(f"processing {genome} ...")
        genome_path = os.path.join(args.input_directory, genome)
        if not os.path.exists(genome_path + ".fai"):
            cmd = ["samtools","faidx",genome_path]
            subprocess.run(cmd)
        lengths = {}
        with open(genome_path + ".fai") as f:
            for line in f:
                seq_id ,length = line.strip().split("\t")[:2]
                length = int(length)
                lengths[seq_id] = length
        logger.info("run minimap2 ...")
        cmd = ["minimap2", "-x", "asm20", genome_path, rep_genome_path]    
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
        handle = io.TextIOWrapper(proc.stdout, encoding="utf-8")
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
            if tname not in intervals:
                intervals[tname] = []
            intervals[tname].append((tstart, tend))
        full_length = 0
        kept_length = 0 
        genome_id = genome[:genome.rfind(".")]
        fasta = pyfaidx.Fasta(genome_path)
        logger.info("extract sequences ...")
        for seq_id in lengths:
            full_length += lengths[seq_id]
            if seq_id not in intervals:
                if lengths[seq_id] < args.min_length:
                    continue
                sequence = str(fasta[seq_id][::])
                fout.write(f">{genome_id}:{seq_id}:0-{lengths[seq_id]}\n")
                p = 0
                while p < len(sequence):
                    fout.write(sequence[p:p+70] + "\n")
                    p += 70
                kept_length += lengths[seq_id]
            else:
                ivs = merge(intervals[seq_id])
                ivs = subtract(ivs, lengths[seq_id])      
                for start, end in ivs:
                    if end - start < args.min_length:
                        continue
                    kept_length += end - start
                    sequence = str(fasta[seq_id][start:end])
                    fout.write(f">{genome_id}:{seq_id}:{start}-{end}\n")
                    p = 0
                    while p < len(sequence):
                        fout.write(sequence[p:p+70] + "\n")
                        p += 70
        logger.info(f"{full_length} {kept_length} {round(kept_length/full_length,4)}")
    fout.close()
    logger.info("All done.")

if __name__ == "__main__":
    main()
