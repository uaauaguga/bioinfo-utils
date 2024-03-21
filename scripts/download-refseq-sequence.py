#!/usr/bin/env python
from Bio import Entrez
import time
import argparse
import os
import logging

def main():
    parser = argparse.ArgumentParser(description='Download genebank sequence')
    parser.add_argument('--query', '-q',help="Genebank id for download",required=True)
    parser.add_argument('--fasta','-f',help="Path to write the downloaded data",required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
    logger = logging.getLogger('Retrieve sequence')

    Entrez.email = "jinyf16@mails.tsinghua.edu.cn"
    fout = open(args.fasta,"w")

    try:
        logger.info("Start retriving {} from ncbi nucleotide in fasta format ...".format(args.query))
        handle = Entrez.efetch(db="nucleotide", id=args.query, rettype="fasta", retmode="text")
        content = handle.read()
        fout.write(content)
        logger.info("Done.")
    except:
        logger.info("Error retriving {}, skip ...".format(args.query))
    fout.close()

if __name__ == "__main__":
    main()
