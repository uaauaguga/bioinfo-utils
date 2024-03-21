#!/usr/bin/env python
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='rename sequence id in fasta file')
    parser.add_argument('--input','-i',type=str, required=True, help='input fasta file')
    parser.add_argument('--output','-o',type=str, required=True, help="output fasta file")
    parser.add_argument('--rename','-lut',type=str, required=True, help="the name mapping file")
    args = parser.parse_args()
    seq_id_lut = {}
    with open(args.rename) as f:
        for line in f:
            seq_id_1, seq_id_2 = line.strip().split("\t")
            seq_id_lut[seq_id_1] = seq_id_2
    fin  = open(args.input)
    if not os.path.exists(args.output): 
        fout = open(args.output,"w")
    else:
        print(args.output + " already exists")
        sys.exit(1)    
    
    for line in fin:
        if line.startswith(">"):
            attrs = line.strip()[1:].split(" ")
            seq_id = attrs[0]
            if seq_id not in seq_id_lut:
                print("{seq_id} not present in {args.rename}, not replace")
            else:
                seq_id = seq_id_lut[seq_id]
            if len(attrs) > 1:
                fout.write(">" + seq_id + " " + " ".join(attrs[1:]) + "\n")
            else:
                fout.write(">" + seq_id + "\n")
        else:
            fout.write(line)

    fin.close()
    fout.close()
            

if __name__ == "__main__":
    main() 
