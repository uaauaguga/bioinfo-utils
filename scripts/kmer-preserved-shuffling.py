#!/usr/bin/env python
import argparse
from ushuffle import shuffle

def main():
    parser = argparse.ArgumentParser(description='shuffle sequence preserve given frequency')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences')
    parser.add_argument('--output','-o',type=str, required=True, help='output sequences')
    parser.add_argument('--kmer','-k',default=2, help="kmer frequency to preserve")
    args = parser.parse_args()
    sequences = {}
    with open(args.input) as f:
        for line in  f:
            if line.startswith(">"):
                header = line
                sequences[header] = ""
            else:
                sequences[header] += line.strip()
    
    with open(args.output,"w") as f:
        for header in sequences:
            sequence = shuffle(sequences[header].encode(),args.kmer).decode()
            f.write(header)
            f.write(sequence + "\n")                    

if __name__ == "__main__":
    main()
