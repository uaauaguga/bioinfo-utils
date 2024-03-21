#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Split Rfam seed alignments, each RNA family as a stockholm format file.')
    parser.add_argument('--input', '-i',required=True,help="Input path")
    parser.add_argument('--outdir','-o',required=True,help="Output path")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    #cope with rfam seed alignment encoding
    with open(args.input,encoding="ISO-8859-1") as f:
        content = ""
        for line in f:
            if line.startswith('#=GF AC'):
                # update accession of the entry
                accession = line.split()[-1]
            elif line.startswith("#=GF ID"):
                name = line.split()[-1]
                accession = accession + "__" + name
            content += line
            if line.startswith("//"):
                fout = open(os.path.join(args.outdir,accession+".stk"),"w")
                fout.write(content)
                fout.close()
                content = ""


if __name__ == "__main__":
    main()

