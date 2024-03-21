#!/usr/bin/env python
import argparse
import gzip
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('chunkify sequences')

def main():
    parser = argparse.ArgumentParser(description='chunkify fasta files ')
    parser.add_argument('--input-directory', '-i', required=True, help="input directory")
    parser.add_argument('--output-directory','-o', required=True, help="output directory")
    parser.add_argument('--chunk-size','-cs', type=int,default=250, help="chunk size to use")
    parser.add_argument('--use-preffix','-p', action="store_true", help="whether use prefix in the output")
    args = parser.parse_args() 
         
    fastas = os.listdir(args.input_directory)
    logger.info(f"{len(fastas)} genomes presented in input directory")

    logger.info(f"use a chunk size of {args.chunk_size} .")
    if not os.path.exists(args.output_directory):
        logger.info(f"output directory {args.output_directory} does not exists, create it")
        os.mkdir(args.output_directory)

    i = 0
    j = 0
    logger.info(f"write to first chunk ...")
    fout = open(os.path.join(args.output_directory,str(i).zfill(4)+".fa"),"w")
    for fasta in sorted(fastas):
        if j >= args.chunk_size:
            i += 1
            j = 0
            logger.info(f"start write to chunk {i} ...")
            fout.close()
            fout = open(os.path.join(args.output_directory,str(i).zfill(4)+".fa"),"w")
        path = os.path.join(args.input_directory, fasta)
        if path.endswith(".gz"):
            f = gzip.open(path)
            prefix = ".".join(fasta.split(".")[:-2])
        else:
            f = open(path)        
            prefix = ".".join(fasta.split(".")[:-1])
        for line in f:
            if path.endswith(".gz"):
                line = line.decode()
            if args.use_preffix and line.startswith(">"):
                line = line.replace(">",f">{prefix}:")
            fout.write(line)
        f.close()
        j += 1    

    try:
        fout.close()
    except:
        pass


if __name__ == "__main__":
    main()
