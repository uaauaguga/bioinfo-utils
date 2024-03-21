#!/usr/bin/env python
import argparse
def main():
    
    parser = argparse.ArgumentParser(description='drop sequence with same sequence id')
    parser.add_argument('--input', '-i', type=str, help='input sequences')
    parser.add_argument('--output','-o', type=str, help='output sequences passed filter')
    args = parser.parse_args() 
    seq_ids = set()
    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                if seq_id in seq_ids:
                    skip = True
                else:
                    skip = False
                    seq_ids.add(seq_id)
            if not skip:
                fout.write(line)
    fout.close()

if __name__ == "__main__":
    main()
            
            
                