#!/usr/bin/env python
import argparse
import numpy as np
import logging
from collections import defaultdict
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("select intervals")



def select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs):
    cached_scores_by_strand = {"+":[],"-":[]}
    cached_iv_by_strand = {"+":[],"-":[]}
    cached_names_by_strand = {"+":[],"-":[]}
    cached_attrs_by_strand = {"+":[],"-":[]}
    for i in range(len(cached_scores)):
        seq_id, start, end, strand = cached_ivs[i]
        cached_scores_by_strand[strand].append(cached_scores[i])
        cached_iv_by_strand[strand].append((seq_id, start, end))
        cached_names_by_strand[strand].append(cached_names[i])
        cached_attrs_by_strand[strand].append(cached_attrs[i])
    score_by_strand, iv_by_strand = {}, {}
    name_by_strand, attrs_by_strand = {}, {}
    level2_cluster_count_by_strand = {}
    for strand in "+-":
        if len(cached_scores_by_strand[strand]) == 0:
            continue
        counter = defaultdict(int)
        for seq_id in cached_names_by_strand[strand]:
            cluster_id_1 = seq_id_to_cluster_id_1[seq_id]
            cluster_id_2 = cluster_id_1_to_cluster_id_2[cluster_id_1]
            counter[cluster_id_2] += 1
        level2_cluster_count_by_strand[strand] = ""
        for cluster_id_2, count in sorted(counter.items(),key=lambda x:-x[1]):
            level2_cluster_count_by_strand[strand] += f"{cluster_id_2}:{count};"
        i = np.argmax(cached_scores_by_strand[strand])
        score_by_strand[strand] = cached_scores_by_strand[strand][i]
        iv_by_strand[strand] = cached_iv_by_strand[strand][i]
        name_by_strand[strand] = cached_names_by_strand[strand][i]
        attrs_by_strand[strand] = cached_attrs_by_strand[strand][i]
    return score_by_strand, iv_by_strand, name_by_strand, attrs_by_strand, level2_cluster_count_by_strand


def main():
    parser = argparse.ArgumentParser(description='pick best interval from overlapped ones')
    parser.add_argument('--input', '-i',required=True,help="input intervals")
    parser.add_argument('--output','-o',required=True,help="output intervals")
    args = parser.parse_args()
    logger.info("load cluster annotation ...")
    global seq_id_to_cluster_id_1
    global cluster_id_1_to_cluster_id_2
    seq_id_to_cluster_id_1 = {}
    cluster_id_1_to_cluster_id_2 = {}
    with open("/apps/home/lulab_jinyunfan/qhsky1/sRNA-finding/output/combined/level2-clusters.txt") as f:
        for line in f:
            seq_id, cluster_id_1, cluster_id_2 = line.strip().split("\t")
            seq_id_to_cluster_id_1[seq_id] = cluster_id_1
            cluster_id_1_to_cluster_id_2[cluster_id_1] = cluster_id_2

    logger.info(f"load intervals from {args.input} ...")
    logger.info(f"picked intervals will be saved to {args.output} .")
    fin = open(args.input)
    seq_id = ""
    last_seq_id, last_start, last_end = "", -1 ,-1
    fout = open(args.output,"w")
    cached_ivs = []
    cached_scores = []
    cached_attrs = []
    cached_names = []
    for line in fin:
        fields = line.strip().split("\t")
        seq_id, start, end, name, score, strand = fields[:6]
        score = float(score)
        attrs = fields[6:]
        start, end = int(start), int(end)
        if ((start > last_end) and (last_seq_id == seq_id)) or (last_seq_id != seq_id):
            # current entry does not overlap with previous one
            # save cached entries
            if len(cached_ivs) > 0:
                # group predictions by strand
                score_by_strand, iv_by_strand, name_by_strand, attrs_by_strand, level2_cluster_count_by_strand = select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs)
                for mstrand in "+-":
                    if mstrand not in score_by_strand:
                        continue
                    mseq_id, mstart, mend  = iv_by_strand[mstrand]
                    mscore = score_by_strand[mstrand]
                    s, e, l = attrs_by_strand[mstrand][:3]
                    s, e, l = int(s), int(e), int(l)
                    if mstrand == "+":
                        mstart = max(mstart - s, 0)
                        mend = mend + (l-e)
                    else:
                        mstart = max(mstart - (l-e), 0)
                        mend = mend + s
                    print(mseq_id, mstart, mend, seq_id_to_cluster_id_1[name_by_strand[mstrand]], round(mscore,3), mstrand, *attrs_by_strand[mstrand], level2_cluster_count_by_strand[mstrand], file=fout, sep="\t")
                # update the cache
                cached_ivs, cached_scores, cached_attrs, cached_names = [], [], [], []
        cached_ivs.append((seq_id, start, end, strand))
        cached_scores.append(score)
        cached_attrs.append(attrs)
        cached_names.append(name)
        last_seq_id, last_start, last_end = seq_id, start, end

    if len(cached_ivs) > 0:
        score_by_strand, iv_by_strand, name_by_strand, attrs_by_strand, level2_cluster_count_by_strand = select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs)
        for mstrand in "+-":
            if mstrand not in score_by_strand:
                continue
            mseq_id, mstart, mend = iv_by_strand[mstrand]
            mscore = score_by_strand[strand]
            s, e, l = attrs_by_strand[mstrand][:3]
            s, e, l = int(s), int(e), int(l)
            if mstrand == "+":
                mstart = max(mstart - s, 0)
                mend = mend + (l-e)
            else:
                mstart = max(mstart - (l-e), 0)
                mend = mend + s
            print(mseq_id, mstart, mend, seq_id_to_cluster_id_1[name_by_strand[mstrand]], round(mscore,3), mstrand, *attrs_by_strand[mstrand], level2_cluster_count_by_strand[mstrand], file=fout, sep="\t")

    fin.close()
    fout.close()

if __name__ == "__main__":
    main()

