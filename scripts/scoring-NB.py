#!/usr/bin/env python
import argparse
import os
import gzip
from itertools import product
import numpy as np
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('NB scoring')

words = {}
rc_lut = {"A":"T","C":"G","G":"C","T":"A"}

def get_kmer_count(sequence, k = 3):
    counts = defaultdict(int)
    for i in range(len(sequence)-k):
        if "N" in sequence[i:i+k]:
            continue
        counts[sequence[i:i+k]] += 1
    return counts

def get_N_count(sequence):
    n_N = 0
    for c in sequence:
        if c == 78:
            n_N += 1
    return n_N

def dust(sequence, k=3):
    if len(sequence) < k:
        return 100
    counts = get_kmer_count(sequence, k=k)
    dust_score = 0
    for count in counts.values():
        dust_score += count*(count-1)
    n_N = get_N_count(sequence)
    if n_N >= 4:
        return 100
    n = len(sequence) - n_N - 2
    ## for polymer of one nucleotide, this score is
    ## n * (n - 1)
    return 100*dust_score/((n - 1)*n)

def scoring(sequence, params, word_size=5):
    score = 0
    n_word = 0
    for i in range(len(sequence)-word_size):
        kmer = sequence[i:i+word_size]
        if kmer in words:
            n_word += 1
            kmer = words[kmer]
            score += np.log(params[kmer])
    return score/n_word        

def main():
    parser = argparse.ArgumentParser(description='Scoring with naive bayes model')
    parser.add_argument('--input-fastq-1', '-if1', type=str, required=True, help='Read 1')
    parser.add_argument('--input-fastq-2', '-if2', type=str, help='Read 2')
    parser.add_argument('--output-fastq-1','-of1',type=str, required=True, help='Output read 1')    
    parser.add_argument('--output-fastq-2','-of2',type=str, help='Output read 2')
    parser.add_argument('--params','-p',type=str, required=True, help='NB parameters')
    args = parser.parse_args()

    
    logger.info("load parameters ...")
    params = {}
    
    with open(args.params) as f:
        header = next(f).strip()
        labels = header.split("\t")[1:]
        for label in labels:
            params[label] = {}
        for line in f:
            fields = line.strip().split("\t")
            kmer = fields[0]
            word_size = len(kmer)
            for label, freq in zip(labels,fields[1:]):
                freq = float(freq)
                params[label][kmer] = freq
    
    for word in list(product(*["ACGT"]*word_size)):
        word = "".join(word)
        word_rc = "".join(rc_lut[c] for c in word[::-1])
        words[word] = min(word, word_rc)
                    
     
    labels = list(params.keys())
    fout = open(args.output_fastq_1,"w")
    if args.output_fastq_2 is not None:
        fout2 = open(args.output_fastq_2,"w")
    #print("seq id",*labels,sep="\t")
    n = 0
    logger.info("start processing ...")
    if args.input_fastq_2 is not None:
        f2 = gzip.open(args.input_fastq_2,"rb")
    with gzip.open(args.input_fastq_1,"rb") as f:
        for header in f:
            n += 1
            if (n%100000) == 0:
                logger.info(f"{round(n/1000000,1)} read processed .")
            header = header.decode().strip()
            sequence = next(f).decode().strip()
            placeholder = next(f).decode().strip()
            quality = next(f).decode().strip()
            score = scoring(sequence, params["viroid"], word_size=word_size)
            sequence2 = ""
            if args.input_fastq_2 is not None:
                header2 = next(f2).decode().strip()
                sequence2 = next(f2).decode().strip()
                placeholder2 = next(f2).decode().strip()
                quality2 = next(f2).decode().strip()
                score2 = scoring(sequence2, params["viroid"], word_size=word_size)
                score = (score+score2)/2
            if score > -6:    
                if (args.input_fastq_2 is not None) and (dust(sequence2) >= 5):
                    continue       
                if dust(sequence) < 5:
                    print(header,file=fout)
                    print(sequence,file=fout)
                    print(placeholder,file=fout)
                    print(quality,file=fout)
                    fout.flush()
                    if args.input_fastq_2 is not None:
                        print(header2,file=fout2)
                        print(sequence2,file=fout2)
                        print(placeholder2,file=fout2)
                        print(quality2,file=fout2)
                        fout2.flush()
    fout.close()
    if args.input_fastq_2 is not None:
        fout2.close()
            

if __name__ == "__main__":
    main()
