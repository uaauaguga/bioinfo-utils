#!/usr/bin/env python
import argparse
import os
from itertools import product
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('fit NB')

words = {}
rc_lut = {"A":"T","C":"G","G":"C","T":"A"}


def get_word_count(path, word_size=5, pseudocount=0.01,max_word=100000000):
    sequence = ""
    counts = {}
    n_word = 0
    with open(path) as f:
        for line in f:
            if n_word > max_word:
                break
            if line.startswith(">"):
                if len(sequence) > 0:
                    for i in range(len(sequence)-word_size):
                        kmer = sequence[i:i+word_size]
                        if kmer in words:
                            n_word += 1
                            kmer = words[kmer]
                            if kmer not in counts:
                                counts[kmer] = 0
                            counts[kmer] += 1
                sequence = ""
            else:
                sequence += line.strip()
    if len(sequence) > 0:
        for i in range(len(sequence)-word_size):
            if n_word > max_word:
                break
            kmer = sequence[i:i+word_size]
            if kmer in words:
                n_word += 1
                kmer = words[kmer]
                if kmer not in counts:
                    counts[kmer] = 0
                counts[kmer] += 1
    for kmer in counts:
        counts[kmer] = (counts[kmer]+pseudocount)/n_word*(1+pseudocount)
    return counts            

def main():
    parser = argparse.ArgumentParser(description='fit naive bayes model')
    parser.add_argument('--input-directory', '-i', type=str, required=True, help='Input fasta files')
    parser.add_argument('--output','-o',type=str, required=True, help='Output word frequencies')
    parser.add_argument('--word-size','-k', type=int , default=5, help="K-mer size to use")
    args = parser.parse_args()
    
    for word in list(product(*["ACGT"]*args.word_size)):
        word = "".join(word)
        word_rc = "".join(rc_lut[c] for c in word[::-1])
        words[word] = min(word, word_rc)
        
    
    frequencies = {}
    for fasta in os.listdir(args.input_directory):
        logger.info(f"processing {fasta} ...")
        label = fasta[:fasta.rfind(".")]
        frequencies[label] = get_word_count(os.path.join(args.input_directory,fasta), word_size=args.word_size)        
    logger.info("saving results ...")
    canonical_words = sorted(list(words.values()))
    labels = sorted(list(frequencies.keys()))
    fout = open(args.output,"w")
    print("word",*labels,sep="\t",file=fout)
    for word in canonical_words:
        freqs = []
        for label in labels:
            freqs.append(frequencies[label].get(word,0))
        print(word,*freqs,sep="\t",file=fout)

    logger.info("all done.")
    

if __name__ == "__main__":
    main()
