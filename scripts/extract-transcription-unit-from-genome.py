#!/usr/bin/env python
import argparse
import numpy as np
import logging
import subprocess
import io
from collections import defaultdict
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("extract TU")

def main():
    parser = argparse.ArgumentParser(description='extract transcription unit from TSS and terminators')
    parser.add_argument('--transcription-start-site', '-tss', type=str, required=True, help='transcription start site')
    parser.add_argument('--transcription-termination-site','-tts', type=str ,required=True, help="transcription termination site")
    parser.add_argument('--transcription-unit','-tu', type=str ,required=True, help="predicted transcription unit")
    parser.add_argument('--min-length','-m', type=int , default=40, help="max length of predicted transcription unit")
    parser.add_argument('--max-length','-M', type=int , default=500, help="max length of predicted transcription unit")
    parser.add_argument('--cutoff', type=float ,  help="cutoff of predicted TSS")
    args = parser.parse_args()

    logger.info("Get closest TSS of each TES ...")
    cmd = ["bedtools","closest","-id","-D","a","-a",args.transcription_termination_site,"-b",args.transcription_start_site,"-s"]
    logger.info("running " + " ".join(cmd))
    proc = subprocess.Popen(cmd,stdout=subprocess.PIPE)
    scores = {}
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        #NC_000913.2     264     314     .       0.645,34.43     -       NC_000913.2     1176    1177    b0002   2.95,3.41,5.27  -       -863
        fields = line.strip().split("\t")
        seq_id = fields[0]
        tss_start, tss_end = int(fields[7]), int(fields[8])
        tts_start, tts_end = int(fields[1]), int(fields[2])
        strand = fields[5]
        if strand == "+":
            start, end = tss_start, tts_end
        else:
            start, end = tts_start, tss_end
        L = end-start
        if (start < 0) or (end < 0):
            continue
        tss_score, tts_score = fields[10], fields[4]
        if (L >= args.min_length) and (L <= args.max_length):
            scores[(seq_id,start,end,strand)] = (tss_score, tts_score)
    logger.info("Get closest TES of each TSS and saving mutually closest pairs ...")
    fout = open(args.transcription_unit,"w")
    cmd = ["bedtools","closest","-iu","-D","a","-b",args.transcription_termination_site,"-a",args.transcription_start_site,"-s"]
    logger.info("running " + " ".join(cmd))
    proc = subprocess.Popen(cmd,stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        fields = line.strip().split("\t")
        tts_start, tts_end = int(fields[7]), int(fields[8])
        tss_start, tss_end = int(fields[1]), int(fields[2])
        strand = fields[5]
        if strand == "+":
            start, end = tss_start, tts_end
        else:
            start, end = tts_start, tss_end
        iv = (seq_id,start,end,strand)
        if iv in scores:
            tss_score, tes_score = scores[iv]
            seq_id, start, end, strand = iv
            tss_score = "." if len(tss_score) == 0 else tss_score
            tts_score = "." if len(tts_score) == 0 else tts_score
            print(seq_id, start, end, tss_score, tes_score, strand,sep="\t", file=fout)
    fout.close()            
          

if __name__ == "__main__":
    main()
