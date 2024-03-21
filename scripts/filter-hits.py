#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("filter hits")

def main():
    parser = argparse.ArgumentParser(description='filter transcript hits')
    parser.add_argument('--input', '-i', type=str, required=True, help='input mmseqs search results')
    parser.add_argument('--output','-o', type=str ,required=True, help="output transcript consensus location")
    parser.add_argument('--min-length','-ml', type=int , default=40, help="minimal alignment length required")
    parser.add_argument('--min-coverage','-mc', type=float , default=0.75, help="minimal alignment coverage required")
    parser.add_argument('--save-bed','-sb', action="store_true", help="whether save results in bed format")
    parser.add_argument('--mutual','-m', action="store_true", help="use mutual coverage instead of query coverage")
    parser.add_argument('--either','-e', action="store_true", help="if specified, only one of min-length/min-coverage needs to be satistified")
    args = parser.parse_args()
    #                     0                                                      1                        2      3      4       5        6      7          8       9        10          11       12       13
    #MGYG-HGUT-00104:GUT_GENOME000545_3:16211-16742:201-461(-)	MGYG-HGUT-00047:GUT_GENOME000200_11	0.869	175	23	0	178	4	136206	136380	1.82E-50	207	260	142274
    fout = open(args.output,"w")
    logger.info("load data ...")
    logger.info(f"hits passed filter will be saved to {args.output} .")
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t") 
            contig_id, tstart, tend = fields[1], int(fields[8]), int(fields[9])
            qstart, qend = int(fields[6]), int(fields[7])
            assert tstart < tend
            query_length, target_length = int(fields[12]), int(fields[13])
            aligned_length = int(fields[3])
            if args.mutual:
                coverage = aligned_length/max(query_length, target_length)
            else:
                coverage = aligned_length/query_length
            if args.either:
                if (coverage < args.min_coverage) and (aligned_length < args.min_length):
                    continue
            else:
                if (coverage < args.min_coverage) or (aligned_length < args.min_length):
                    continue
            if args.save_bed:
                if qstart < qend:
                    strand = "+"
                else:
                    strand = "-"
                    qstart, qend = qend, qstart
                qstart -= 1
                coverage = round(coverage,3)
                identity = round(float(fields[2]),3)
                print(contig_id,tstart-1,tend,fields[0],f"{coverage*identity}",strand, qstart, qend, query_length, identity, sep="\t",file=fout)
            else:
                fout.write(line)
    fout.close()
    logger.info("all done .")
    
    
            
            


if __name__ == "__main__":
    main()
