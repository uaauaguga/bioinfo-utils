#!/usr/bin/env python
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Split Rfam seed alignments, each RNA family as a stockholm format file.')
    parser.add_argument('--input', '-i',required=True,help="Input path")
    parser.add_argument('--outdir','-o',required=True,help="Output path")
    parser.add_argument('--use-name','-un', action = 'store_true',help="Use name instaed of accession")
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    
    with open(args.input) as f:
        content = ""
        for line in f:
            if line.startswith('NAME'):
                # update accession of the entry
                name = line.split()[-1]
            elif line.startswith("ACC"):
                accession = line.split()[-1]
                accession = accession + "__" + name
            content += line
            if line.startswith("//"):
                if not content.strip().startswith("INFERNAL"):
                    print(name)
                    if not args.use_name:
                        path = os.path.join(args.outdir,accession+".cm")
                    else:
                        path = os.path.join(args.outdir,name+".cm")
                    fout = open(path,"w")
                    fout.write(content0)
                    fout.write(content)
                    fout.close()
                else:
                    content0 = content
                content = ""


if __name__ == "__main__":
    main()

