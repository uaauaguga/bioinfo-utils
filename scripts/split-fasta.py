#!/usr/bin/env python
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='split fasta file to multiple one')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='directory to save chunked fasta')
    parser.add_argument('--chunk-size','-cs', type=int,default=1000, help="chunk size to use")
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    i = 0
    n = 0
    fout = open(os.path.join(args.output_directory,str(i).zfill(4) + ".fa"),"w")
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                # new sequence, check whether open a new file
                n += 1
                if n >= args.chunk_size:
                    fout.close()
                    i += 1
                    fout = open(os.path.join(args.output_directory,str(i).zfill(4) + ".fa"),"w")
                    n = 0
            fout.write(line)
    try:
        fout.close()
    except:
        pass

if __name__ == "__main__":
    main()
