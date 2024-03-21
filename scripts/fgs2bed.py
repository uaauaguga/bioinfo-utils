#!/usr/bin/env python
import argparse
import re
def main():
    parser = argparse.ArgumentParser(description='convert fraggenescan to bed format')
    parser.add_argument('--input', '-i', type=str, required=True, help='input fraggenescan results')
    parser.add_argument('--output','-o', type=str ,required=True, help="output in bed format")
    args = parser.parse_args()
    # >k119_9373:1630-2842
    # 128	541	-	2	1.375697	I:	D:
    # 1031	1210	-	2	1.402302	I:	D: 
    fout = open(args.output,"w")
    with open(args.input) as fin:
        for line in fin:
            if line.startswith(">"):
                seq_id = line.strip()[1:]
            else:
                fields = re.split(r"\s+",line.strip())        
                start, end, strand, score = int(fields[0])-1, int(fields[1]), fields[2], float(fields[4])
                print(seq_id, start, end,".",score, strand, file=fout,sep="\t")


if __name__ == "__main__":
    main()
