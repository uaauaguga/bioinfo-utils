#!/usr/bin/env python
import sys
import subprocess
import argparse
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('prepare shuffled sequences')

def main():
    parser = argparse.ArgumentParser(description='preapre shuffled sequences')
    parser.add_argument('--input',  '-i', type=str, required=True, help="Input sequence in fasta format")
    parser.add_argument('--output-prefix', '-op', type=str, required=True, help="Output shuffled sequence in fasta format")
    parser.add_argument('--number', '-n', type=int, default=50, help="Number of shuffling")
    args = parser.parse_args()        

    logger.info("rename input sequence ...")
    frename = open(args.input + ".renamed","w")
    i = 0
    seq_ids = []
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split(" ")[0]
                seq_ids.append(seq_id)
                line = f">{i}\n" 
                i += 1
            frename.write(line)
    frename.close()
                
    
    logger.info("align input sequence ...")
    falignment = open(args.input + ".aligned.aln","w")
    cmd = ["mafft","--auto", "--clustalout", args.input + ".renamed"]
    proc = subprocess.run(cmd, stdout=falignment, stderr = subprocess.DEVNULL)
    falignment.close()     

    shuffled_path = args.input + ".aligned.shuffled.fa"
    logger.info("shuffling ...")
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/tools/SISSIz/src/SISSIz", "-n", 
           str(args.number), "--fasta", "-s", "-o", shuffled_path, args.input + ".aligned.aln"]
    print(" ".join(cmd))
    subprocess.run(cmd)

    logger.info("prepare shuffled sequence ...")
    chunks = open(shuffled_path).read().split("//\n")
    i = 0

    logger.info("summarize results ...")
    for chunk in chunks:
        if len(chunk) == 0:
            continue
        sequences = {}
        for line in chunk.split("\n"):
            if len(line) == 0:
                continue
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip().replace("-","")
        fout = open(args.output_prefix + f".{i}.fa","w")
        for seq_id in sorted(sequences.keys()):
            print(f">shuffled-{seq_ids[int(seq_id)]}",file=fout)
            print(sequences[seq_id],file=fout)
        fout.close()
        i += 1
    fout.close()

    logger.info("clean up ...")
    os.remove(args.input + ".aligned.aln")
    os.remove(shuffled_path)
    os.remove(args.input + ".renamed")
    logger.info("all done.")

if __name__ == "__main__":
    main()
