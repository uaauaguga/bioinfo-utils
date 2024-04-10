#!/usr/bin/env python
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='split fasta file to multiple one')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='directory to save chunked fasta')
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    last_genome_id = ""
    fout = None
    with open(args.input) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split(" ")[0]
                genome_id = seq_id[:seq_id.find(":")]
                line = ">" + line[line.find(":")+1:]
                if last_genome_id != genome_id:
                    if fout:
                        fout.close()
                    fout = open(os.path.join(args.output_directory,f"{genome_id}.fna"),"w")    
                last_genome_id = genome_id
            fout.write(line)
    try:
        fout.close()
    except:
        pass

if __name__ == "__main__":
    main()
