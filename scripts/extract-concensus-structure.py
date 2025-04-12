#!/usr/bin/env python
import argparse
import logging
from collections import defaultdict
import os
from Bio import AlignIO
import numpy as np
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract consensus structure')

unpaired_chars = "_-.:,~"
left_chars =  "<([{ABCDE"
right_pairs = ">)]}abcde"
forward_lookup = dict(zip(left_chars,right_pairs))
backward_lookup = dict(zip(right_pairs,left_chars))
base_lookup = {"A":"U","C":"G","G":"CU","U":"AG","N":""}

def load_consensus_pairs(consensus_structure):
    stacks = defaultdict(list)
    pseudoknot_consensus_pairs = defaultdict(list)
    pseudoknot_free_consensus_pairs = []
    for p, c in enumerate(consensus_structure):
        if c in unpaired_chars:
            continue
        if c in forward_lookup:
            stacks[c].append(p)
        if c in backward_lookup:
            pl, pr = stacks[backward_lookup[c]].pop(), p
            if c not in "abcd":
                pseudoknot_free_consensus_pairs.append((pl, pr))
            else:
                pseudoknot_consensus_pairs[c.upper()].append((pl, pr))
    for c in stacks:
        if len(stacks[c]) != 0:
            print("unpaired base detected, the structure is not valid")
            return None,None
    return sorted(pseudoknot_free_consensus_pairs), pseudoknot_consensus_pairs

def main():
    parser = argparse.ArgumentParser(description='extract consensus structure for stockhol format file')
    parser.add_argument('--input', '-i', type=str, required=True, help='input MSA with structure annotation')
    parser.add_argument('--output','-o', type=str , required = True, help="where to save output")
    parser.add_argument('--stat','-s', type=str , required=True, help="where to save statistics")
    args = parser.parse_args()
    logger.info("load MSA ...")
    msa = AlignIO.read(args.input,format="stockholm")
    logger.info("load consensus pairing ...")
    consensus_structure = msa.column_annotations['secondary_structure']
    pseudoknot_free_consensus_pairs, pseudoknot_consensus_pairs = load_consensus_pairs(consensus_structure)
    fout  = open(args.output,"w")
    fstat = open(args.stat,"w")
    logger.info("extract structure of each sequence ...")
    for record in msa:
        seq_id, sequence = record.id, str(record.seq)
        ungapped_sequence = sequence.replace("-","")
        mask = np.array(list(sequence)) != "-"
        cp2sp = mask.cumsum() - 1
        n_gapped_pairs, n_mismatched_pairs = 0, 0
        n_pseudoknot_pairs, n_nested_pairs = 0, 0
        dbn = ["." for i in range(len(ungapped_sequence))]
        consensus_pairs = defaultdict()
        consensus_pairs["free"] = pseudoknot_free_consensus_pairs
        consensus_pairs.update(pseudoknot_consensus_pairs)
        for t in consensus_pairs:
            for lp, rp in consensus_pairs[t]:
                slp, srp = cp2sp[lp], cp2sp[rp]
                if (sequence[lp] == "-") or (sequence[rp] == "-"):
                    n_gapped_pairs += 1
                elif sequence[rp] in base_lookup.get(sequence[lp],"N"):
                    assert ungapped_sequence[srp] in base_lookup.get(ungapped_sequence[slp])
                    if t == "free":
                        dbn[slp], dbn[srp] = "(", ")"
                        n_nested_pairs += 1
                    else:
                        dbn[slp], dbn[srp] = t, t.lower()
                        n_pseudoknot_pairs += 1
                else:
                    n_mismatched_pairs += 1
        print(seq_id,n_nested_pairs,n_pseudoknot_pairs, n_mismatched_pairs, n_gapped_pairs, 
              round(2*(n_nested_pairs+n_pseudoknot_pairs)/len(ungapped_sequence),3), sep="\t", file=fstat)
        print(f">{seq_id}",file=fout)
        print(ungapped_sequence,file=fout)
        dbn = "".join(dbn)
        print(dbn,file=fout)
    fout.close()
    fstat.close()

if __name__ == "__main__":
    main()

    
            
      
       
