#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='MCL output to clustering table')
    parser.add_argument('--input', '-i', type=str, required=True, help='input MCL output')
    parser.add_argument('--output','-o', type=str ,required=True, help="output clustering table")
    args = parser.parse_args() 
    fout = open(args.output,"w")
    with open(args.input) as fin:
        for i,line in enumerate(fin):
            fields = line.strip().split("\t")
            cluster_id = str(i).zfill(10)
            for seq_id in fields:
                print(seq_id, cluster_id, len(fields), sep="\t",file=fout)
    fout.close()


if __name__ == "__main__":
    main()
