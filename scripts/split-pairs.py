#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='split sequence')
    parser.add_argument('--input', '-i', type=str, required=True, help='predicted  interaction')
    parser.add_argument('--first-rna',  '-fr', type=str, required=True, help='first RNA sequences')
    parser.add_argument('--second-rna', '-sr', type=str, required=True, help='second RNA sequences')
    args = parser.parse_args()

    fin = open(args.input)
    f1 = open(args.first_rna,"w")
    f2 = open(args.second_rna,"w")

    for header in fin:
        sequence = next(fin)
        if sequence.count("NNN") != 1:
            continue
        sequence_1, sequence_2 = sequence.strip().split("NNN")
        print(header.strip(),file=f1)
        print(sequence_1,file=f1)
        print(header.strip(),file=f2)
        print(sequence_2,file=f2)
    f1.close()
    f2.close()
    fin.close()


if __name__ == "__main__":
    main()
