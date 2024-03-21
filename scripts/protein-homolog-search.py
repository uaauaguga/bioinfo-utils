#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("pick hits")
import subprocess
import os
from collections import defaultdict
from pyfaidx import Fasta
from Bio.Seq import Seq

def main():
    parser = argparse.ArgumentParser(description='protein homolog search')
    parser.add_argument('--query', '-q',required=True,help="input query proteins")
    parser.add_argument('--genome-directory','-gd',required=True,help="target genomes")
    parser.add_argument('--cds-directory','-cd',required=True,help="CDS directory")
    parser.add_argument('--output-directory','-od',required=True,help="output directory")
    parser.add_argument('--coverage','-c',default=0.8,help="alignment coverage required")
    parser.add_argument('--threads','-t',default=8,help="threads for mmseqs search")
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    genome_fasta = os.path.join(args.output_directory,"genomes.fa")
    if not os.path.exists(genome_fasta):
        logger.info("prepare genome sequences ...")
        fout = open(genome_fasta,"w")
        for fasta in os.listdir(args.genome_directory):
            if not (fasta.endswith(".fa") or fasta.endswith(".fna")):
                continue
            path = os.path.join(args.genome_directory,fasta)
            genome_id = fasta[:fasta.rfind(".")]
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        line = line[:1] + genome_id + ":" + line[1:]
                    fout.write(line)
        fout.close()


    coordinate_lookup = {}
    logger.info("extract protein sequences ...")
    fasta = Fasta(genome_fasta)
    protein_fasta = os.path.join(args.output_directory,"proteins.faa")
    fout = open(protein_fasta,"w")
    for fa in os.listdir(args.genome_directory):
        genome_id = fa[:fa.rfind(".")]
        bed = os.path.join(args.cds_directory,genome_id+".bed")
        with open(bed) as f:
            for line in f:
                fields = line.strip().split("\t")
                seq_id, start, end, protein_id, _, strand = fields[:6]
                start, end = int(start), int(end)
                seq_id = genome_id + ":" + seq_id
                protein_id = genome_id + ":" + protein_id
                start, end = int(start), int(end)
                sequence = fasta[seq_id][start:end]
                if strand == "-":
                    sequence = sequence.reverse.complement
                protein = str(Seq(str(sequence)).translate())
                coordinate_lookup[protein_id] = (seq_id,start,end,strand)
                print(f">{protein_id}",file=fout)
                print(protein,file=fout)
    fout.close()


    tdb = os.path.join(args.output_directory,"proteome")
    if not os.path.exists(tdb + ".dbtype"):
        logger.info("build proteome mmseqs database ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/mmseqs-env/bin/mmseqs", "createdb", "--dbtype", "1", protein_fasta, tdb ]
        subprocess.run(cmd)
    else:
        logger.info("proteome db exists .")

    qdb = os.path.join(args.output_directory,"proteins")
    if not os.path.exists(qdb + ".dbtype"):
        logger.info("build protein mmseqs database ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/mmseqs-env/bin/mmseqs", "createdb", "--dbtype", "1", args.query, qdb ]
        subprocess.run(cmd)
    else:
        logger.info("protein db exists .")

    hdb = os.path.join(args.output_directory,"hits")
    if not os.path.exists(hdb+".dbtype"):
        logger.info("homolog search ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/mmseqs-env/bin/mmseqs", "search", "--search-type", "1","-s", "7.5", "--max-seqs", "5000", qdb, tdb, hdb, "tmp" ]
        subprocess.run(cmd)
    else:
        logger.info("hit db exists .")


    tsv = os.path.join(args.output_directory,"hits.tsv")
    if not os.path.exists(tsv):
        logger.info("reformat homolog search hits ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/mmseqs-env/bin/mmseqs", "convertalis", "--format-mode", "2", qdb, tdb, hdb, tsv ]
        subprocess.run(cmd)
    else:
        logger.info("hit file exists .")
            
    logger.info("extract best hits by genome ...")
    bitscore_by_genome = defaultdict(int)
    best_hits = {}
    with open(tsv) as f:
        for line in f:
            fields = line.strip().split("\t")
            hit_id, hstart, hend = fields[1], int(fields[8]), int(fields[9])
            query_id, qstart, qend = fields[0], int(fields[6]), int(fields[7])
            query_length, hit_length = int(fields[12]), int(fields[13])
            aligned_length = int(fields[3])
            coverage = aligned_length/query_length
            if coverage < args.coverage:
                continue     
            genome_id = hit_id[:hit_id.find(":")]
            bitscore = int(fields[11]) 
            if bitscore > bitscore_by_genome[(query_id, genome_id)]:
                seq_id,start,end,strand = coordinate_lookup[hit_id] 
                if strand == "+":
                    start, end = start - 200, start + 100
                else:
                    start, end = end - 100, end + 200
                if (start < 0) or (end > len(fasta[seq_id])):
                    continue
                best_hits[(query_id,genome_id)] = (seq_id, start, end, strand)
                bitscore_by_genome[(query_id, genome_id)] = bitscore
   
    fout = open(os.path.join(args.output_directory,"leaders.fa"),"w")
    logger.info("save homolog sequences ...")    
    for query_id,genome_id in sorted(list(best_hits.keys())):
        contig_id, start, end, strand = best_hits[(query_id,genome_id)]
        sequence = fasta[contig_id][start:end]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)
        print(f">{query_id}--{contig_id}:{start}-{end}({strand})",file=fout)
        print(sequence,file=fout)
    fout.close()  

              

if __name__ == "__main__":
    main()
