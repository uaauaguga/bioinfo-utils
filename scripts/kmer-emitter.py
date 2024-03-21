#!/usr/bin/env python
from collections import defaultdict
import argparse
import logging
import numpy as np
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("simulate sequence")

def main():
    parser = argparse.ArgumentParser(description='simulate sequence given kmer frequency')
    parser.add_argument('--frequency',  '-f', type=str, required=True, help='input fasta')
    parser.add_argument('--output', '-o', type=str, required=True, help='output frequency table')
    parser.add_argument('--length', '-l', type=int, default=1000, help='length of simulated sequence')
    parser.add_argument('--number', '-n', type=int, default=1000, help='number of sequence to simulate')
    parser.add_argument('--prefix', '-p', help='prefix to prepend to the sequence id line')
    parser.add_argument('--seed', '-s', type=int, default=666, help='random seed to use')
    args = parser.parse_args()

    np.random.seed(args.seed)

    logger.info("load k-mer frequency ...")
    frequencies = {}

    word_size = None
    scaler = 0
    with open(args.frequency) as f:
        for line in f:
            if line.startswith("#"):
                continue
            word, frequency = line.strip().split("\t")
            frequency = float(frequency)
            frequencies[word]  = frequency
            scaler += frequency
            if word_size is None:
                word_size = len(word)
            else:
                assert word_size == len(word)
    words = list(sorted(frequencies.keys()))
    
    logger.info("get next character distribution ...")
    # given k-1 mer, look for the next character

    conditioned_distribution0 = defaultdict(dict)
    for word in frequencies:
        conditioned_distribution0[word[:-1]][word[-1]] = frequencies[word]

    conditioned_distribution = {}
    for condition in conditioned_distribution0:
        distribution = []
        d = conditioned_distribution0[condition]
        for c in "ACGT":
            if c in d:
                distribution.append(d[c])
            else:
                distribution.append(0)
        conditioned_distribution[condition] = np.array(distribution)/np.array(distribution).sum()
    del conditioned_distribution0
    frequencies = [ frequencies[word] for word in words ]
    frequencies = np.array(frequencies)
    frequencies = frequencies/frequencies.sum()
    logger.info("simulate sequences ...")

    logger.info(f"results are saving to {args.output} ...")
    fout = open(args.output,"w")

    prefix = args.prefix  + ":" if args.prefix is not None else ""

    width = max(int(np.log10(args.number) + 2),4)
    characters = list("ACGT")
    for i in range(args.number):
        if (i+1)%500 == 0:
            logger.info(f"{round(i/1000,2)} K sequence simulated .")
        sequence = np.random.choice(words,p=frequencies)
        condition = sequence[1:]
        index = str(i).zfill(width) 
        header = f">{prefix}{index}:{sequence}"
        fout.write(header + "\n")
        for j in range(args.length-word_size):
            c = np.random.choice(characters, p=conditioned_distribution[condition])
            sequence += c
            condition = condition[1:] + c
        p = 0
        while p < len(sequence):
            fout.write(sequence[p:p+70] + "\n")
            p += 70
    fout.close()

    logger.info("all done .")


if __name__ == "__main__":
    main()



            
            
