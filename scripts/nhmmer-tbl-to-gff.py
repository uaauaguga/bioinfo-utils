#!/usr/bin/env python
import argparse
import re
import numpy as np

def attr_formatter(attrs):
    s = ""
    for k in attrs:
        v = attrs[k]
        if v == "-":
            if k == "ID" and "RNA_name" in attrs:
                v = attrs["RNA_name"]
            else:
                v = "None"
        s += f"{k}={v};"
    return s

def main():
    parser = argparse.ArgumentParser(description='convert nhmmer output to gff format')
    parser.add_argument('--input',  '-i', required = True, help="input nhmmer hits table")
    parser.add_argument('--output', '-o', required = True, help="output hits in gff format")
    args = parser.parse_args()


    fout = open(args.output,"w")
    with open(args.input) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or len(line) == 0:
                continue
            fields = re.split("\s+",line) 
            # 0.  target seq name
            # 1.  target accession
            # 2.  query profile name 
            # 3.  query profile acession
            # 4.  hmmfrom
            # 5.  hmm to
            # 6.  alifrom
            # 7.  ali to  
            # 8.  envfrom
            # 9.  env to
            # 10. sq len
            # 11. strand
            # 12. evalue
            # 13. score 
            # 14. bias
            # 15. description of target

            qseq_id, tseq_id = fields[2], fields[0]
            qseq_acc, tseq_acc = fields[3], fields[1]
            qstart, qend = int(fields[4]), int(fields[5])
            tstart, tend = int(fields[6]), int(fields[7])
            
            assert qstart <  qend

            strand = "+" if (tstart < tend) else "-"
            assert strand == fields[11]
            e_value = fields[12]
            score = fields[13]
            coverage = round(100*(qend - qstart)/float(fields[10]),5)
            attrs = {"ID":qseq_acc,
                     "RNA_name":qseq_id,
                     "e_value":e_value,
                     "coverage":coverage}
            print(tseq_id, 
                  "nhmmer","nucleotide", 
                  tstart,tend,score,
                  strand, "0",
                  attr_formatter(attrs), file=fout, sep="\t") 
    fout.close()


if __name__ == "__main__":
    main()
            
