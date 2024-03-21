#!/usr/bin/env python
import logging
import argparse
import os
from itertools import chain
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('categorization')

def parseAttr(s):
    info = {}
    s = s.strip()
    for data in s.split(";"):
        data = data.strip()
        if len(data) == 0:
            continue
        if "=" in data:
            key,value = data.split("=")
        else:
            key,value = data.split(" ")
        key = key.replace('"','')
        value = value.replace('"','')
        info[key] = value
    return info

def subtract(iv, ivs):
    # subtract intervals from a interval
    start, end = iv
    ivs0 = sorted(ivs)
    ivs = []
    for s, e in ivs0:
        if e < start:
            continue
        if s < start:
            s = start
        if s > end:
            break
        if e > end:
            e = end
        if s == e:
            continue
        ivs.append((s,e))
    ivs = list(chain.from_iterable([(start,)]+ivs+[(end,)]))
    it = iter(ivs)
    subtracted_ivs = []
    for s in it:
        e = next(it)
        if s != e:
            subtracted_ivs.append((s,e))
    return subtracted_ivs

def merge(ivs):
    # merge overlapping intervals
    s, e = -1, -1
    merged_ivs = []
    for start, end in sorted(ivs) + [(1e10000,1e10000)]:
        if start > e:
            # cannot merge, save the entry start a new interval
            if e > 0:
                merged_ivs.append((s,e))
            s = start
        e = max(e,end)
    return merged_ivs

def parse_last_gene():
    global exons_by_tx,cdss_by_tx,introns_by_tx,five_prime_utrs_by_tx,three_prime_utrs_by_tx
    if len(exons_by_tx) > 0:
        exons_by_gene = merge([iv for ivs in exons_by_tx.values() for iv in ivs])
        for ebgs, ebge in exons_by_gene:
            print(chrom, ebgs, ebge, f"{gene_type}|{gene_id}|{gene_name}", ".", strand,sep="\t",file=f_gene_exon)
        introns_by_gene = subtract((gene_start, gene_end),exons_by_gene)
        for ibgs, ibge in introns_by_gene:
            print(chrom, ibgs, ibge, f"{gene_type}|{gene_id}|{gene_name}", ".", strand,sep="\t",file=f_gene_intron)
    else:
        # the last gene contains no transcripts?
        pass
    if len(cdss_by_tx) > 0:
        cdss_by_gene = merge([iv for ivs in cdss_by_tx.values() for iv in ivs])
        for cs, ce in cdss_by_gene:
            print(chrom, cs, ce, f"{gene_type}|{gene_id}|{gene_name}", ".", strand,sep="\t",file=f_gene_cds) 
        if len(five_prime_utrs_by_tx) > 0:
            merged_five_prime_utr = sorted(merge(list(five_prime_utrs_by_tx.values())))
            merged_five_prime_utr = merged_five_prime_utr[0][0],merged_five_prime_utr[-1][1]
            five_prime_utrs_by_gene = subtract(merged_five_prime_utr,[(cdss_by_gene[0][0],cdss_by_gene[-1][1])])
            if len(five_prime_utrs_by_gene) > 0:
                if strand == "+":
                    print(chrom, five_prime_utrs_by_gene[0][0], five_prime_utrs_by_gene[0][1], 
                      f"{gene_type}|{gene_id}|{gene_name}", ".",strand,sep="\t",file=f_gene_5pUTR)
                else:
                    print(chrom, five_prime_utrs_by_gene[-1][0], five_prime_utrs_by_gene[-1][1], 
                      f"{gene_type}|{gene_id}|{gene_name}", ".",strand,sep="\t",file=f_gene_5pUTR)
        if len(three_prime_utrs_by_tx) > 0:
            merged_three_prime_utr = sorted(merge(list(three_prime_utrs_by_tx.values())))
            merged_three_prime_utr = merged_three_prime_utr[0][0],merged_three_prime_utr[-1][1]                    
            three_prime_utr_by_gene = subtract(merged_three_prime_utr,[(cdss_by_gene[0][0],cdss_by_gene[-1][1])])
            if len(three_prime_utr_by_gene) > 0:
                if strand == "+":
                    print(chrom, three_prime_utr_by_gene[-1][0], three_prime_utr_by_gene[-1][1], 
                      f"{gene_type}|{gene_id}|{gene_name}", ".",strand,sep="\t",file=f_gene_3pUTR)     
                else:
                    print(chrom, three_prime_utr_by_gene[0][0], three_prime_utr_by_gene[0][1], 
                      f"{gene_type}|{gene_id}|{gene_name}", ".",strand,sep="\t",file=f_gene_3pUTR)    
    exons_by_tx,cdss_by_tx,introns_by_tx,five_prime_utrs_by_tx,three_prime_utrs_by_tx = {}, {}, {}, {}, {}

def parse_last_tx():
    global exons, exons_by_tx, three_prime_utrs_by_tx, five_prime_utrs_by_tx, cdss
    if len(exons) > 0:
        #print(exons)
        # summarize last transcript
        exons_by_tx[tx_id] = exons
        for es, ee in exons:
            print(chrom, es, ee, f"{gene_type}|{gene_id}|{gene_name}|{tx_id}", ".", strand,sep="\t",file=f_tx_exon) 
        introns_by_tx[tx_id] = subtract((tx_start, tx_end), exons)
        for ins, ine in  introns_by_tx[tx_id]:
            print(chrom, ins, ine, f"{gene_type}|{gene_id}|{gene_name}|{tx_id}", ".", strand,sep="\t",file=f_tx_intron) 
        five_prime_utr, three_prime_utr = None, None
        #print(cdss)
        if len(cdss) > 0:
            for cs, ce in cdss:
                print(chrom, cs, ce, f"{gene_type}|{gene_id}|{gene_name}|{tx_id}", ".", strand,sep="\t",file=f_tx_cds) 
            cdss_by_tx[tx_id] = cdss
            if tx_start < cds_start:
                five_prime_utr = (tx_start, cds_start)
            if cds_end < tx_end:
                three_prime_utr = (cds_end,tx_end)
            if strand == "-":
                five_prime_utr, three_prime_utr = three_prime_utr, five_prime_utr                                        
            if five_prime_utr is not None:
                five_prime_utrs_by_tx[tx_id] = five_prime_utr
                print(chrom, *five_prime_utr, f"{gene_type}|{gene_id}|{gene_name}|{tx_id}", ".", strand,sep="\t",file=f_tx_5pUTR)  
            if three_prime_utr is not None:
                three_prime_utrs_by_tx[tx_id] = three_prime_utr
                print(chrom, *three_prime_utr, f"{gene_type}|{gene_id}|{gene_name}|{tx_id}", ".", strand,sep="\t",file=f_tx_3pUTR) 
    else:
        # last transcript contains no exon?
        pass    
    exons, cdss = [], []


exons, introns, five_prime_utr, three_prime_utr = [], [], None, None
exons_by_tx, cdss_by_tx, introns_by_tx, five_prime_utrs_by_tx, three_prime_utrs_by_tx = {}, {}, {}, {}, {}
chrom, gene_id, gene_name, gene_type, strand = "", "", "", "", ""

def main():
    parser = argparse.ArgumentParser(description='categorize gene by annotation, written in pure python with no dependency')
    parser.add_argument('--gff','-g',type=str, required=True, help="input annotation in gff format")
    parser.add_argument('--output-directory','-od',type=str, required=True, help="categorized genomic intervals")
    parser.add_argument('--biotype','-b',type=str, default="rename.txt", help="how to summarize the biotypes")
    args = parser.parse_args()


    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    global f_gene_5pUTR, f_gene_cds, f_gene_3pUTR, f_gene_intron, f_gene_exon
    global f_tx_5pUTR, f_tx_cds, f_tx_3pUTR, f_tx_intron, f_tx_exon
    global tx_id, tx_start, tx_end, cds_start, cds_end, gene_start, gene_end, gene_id, gene_name, gene_type, chrom, strand
    global exons, exons_by_tx, three_prime_utrs_by_tx, five_prime_utrs_by_tx, cdss,cdss_by_tx,introns_by_tx

    # output files 
    # exon by transcripts
    f_tx_exon = open(os.path.join(args.output_directory,"tx.exon.bed"),"w")
    # CDS by transcripts
    f_tx_cds = open(os.path.join(args.output_directory,"tx.CDS.bed"),"w")
    # UTR by transcripts
    f_tx_5pUTR = open(os.path.join(args.output_directory,"tx.5pUTR.bed"),"w")
    f_tx_3pUTR = open(os.path.join(args.output_directory,"tx.3pUTR.bed"),"w")
    # intron by transcripts
    f_tx_intron = open(os.path.join(args.output_directory,"tx.intron.bed"),"w")

    # exon by genes: exon in at least 1 isoform        
    f_gene_exon = open(os.path.join(args.output_directory,"gene.exon.bed"),"w")
    # CDS by genes: at least is CDS in 1 isoform 
    f_gene_cds = open(os.path.join(args.output_directory,"gene.CDS.bed"),"w")
    # UTR by gene: Not CDS in any isoform
    f_gene_5pUTR = open(os.path.join(args.output_directory,"gene.5pUTR.bed"),"w")
    f_gene_3pUTR = open(os.path.join(args.output_directory,"gene.3pUTR.bed"),"w")
    # intron by genes: not exon in all isoform
    f_gene_intron = open(os.path.join(args.output_directory,"gene.intron.bed"),"w")
        
    logger.info("load annotations ...")
    with open(args.gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            feature = fields[2]
            attrs = parseAttr(fields[8])
            # processing annotations, one gene at a time
            if feature == "gene":
                # save gene model of cahched transcripts, if any
                parse_last_tx()
                parse_last_gene()
                # encounter a new gene
                gene_id, gene_name, gene_type = attrs["gene_id"],attrs["gene_name"],attrs["gene_type"]
                #gene_type = rename_lookup[gene_type]
                # always use 0 base coordinate
                gene_start, gene_end = int(fields[3])-1, int(fields[4])
                strand = fields[6]
            
            elif feature == "transcript":
                # encounter a new transcript
                parse_last_tx()
                assert attrs["gene_id"] == gene_id
                tx_id = attrs["transcript_id"]
                tx_start, tx_end = int(fields[3])-1, int(fields[4])
                
                
            elif feature == "exon":
                assert attrs["gene_id"] == gene_id
                assert attrs["transcript_id"] == tx_id
                exon_start, exon_end = int(fields[3])-1, int(fields[4])
                exons.append((exon_start, exon_end))
                
            elif feature == "CDS":
                assert attrs["gene_id"] == gene_id
                assert attrs["transcript_id"] == tx_id
                cds_start, cds_end = int(fields[3])-1, int(fields[4])
                cdss.append((cds_start, cds_end))

    parse_last_tx()
    parse_last_gene()
    f_gene_5pUTR.close()
    f_gene_cds.close()
    f_gene_3pUTR.close()
    f_gene_intron.close()
    f_gene_exon.close()
    f_tx_5pUTR.close()
    f_tx_cds.close()
    f_tx_3pUTR.close()
    f_tx_intron.close()
    f_tx_exon.close()

if __name__ == "__main__":
    main()
