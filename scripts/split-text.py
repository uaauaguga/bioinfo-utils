#!/usr/bin/env python
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='split text files into multiple ones')
    parser.add_argument('--input', '-i', type=str, required=True, help='input file in bed format')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='directory to save bed files')
    parser.add_argument('--chunk-size','-cs', type=int,default=10000, help="chunk size to use")
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    i = 0
    n = 0
    # the width should be sufficiently large
    fout = open(os.path.join(args.output_directory,str(i).zfill(6) + ".txt"),"w")
    with open(args.input) as f:
        for line in f:
            n += 1
            fout.write(line)
            if n > args.chunk_size:
                fout.close()
                i += 1
                fout = open(os.path.join(args.output_directory,str(i).zfill(4) + ".txt"),"w")
                n = 0
    try:
        fout.close()
    except:
        pass

if __name__ == "__main__":
    main()
