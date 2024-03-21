#!/usr/bin/env python
import argparse
import os
import numpy as np
import subprocess
import logging
import re
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('cmfinder search')


rerun_all = True

def load_fasta(path):
    sequences = {}
    attrs = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                fields = line[1:].strip().split(" ")
                seq_id = fields[0]
                attr = fields[1:]
                sequences[seq_id] = ""
                attrs[seq_id] = attr
            else:
                sequences[seq_id] += line.strip()
    return sequences, attrs

cmfinder_exe = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/cmfinder-env/bin/cmfinder04.pl"
cmbuild_exe = "/apps/home/lulab_jinyunfan/qhsky1/terminator-prediction/terminator-annotation/tools/infernal-1.0.2/src/cmbuild"
cmsearch_exe = "/apps/home/lulab_jinyunfan/qhsky1/terminator-prediction/terminator-annotation/tools/infernal-1.0.2/src/cmsearch" 
rscape_exe = "/apps/home/lulab_jinyunfan/qhsky1/tools/rscape_v2.0.0.j/bin/R-scape"

def cmfinder():
    fasta_path = os.path.join(args.output_directory, name + ".90.sampled.fa")
    log_path = os.path.join(args.output_directory, name + ".90.sampled.cmfinder.log")
    flog = open(log_path,"wb")
    cwd = args.output_directory
    fasta = os.path.basename(fasta_path)
    indir = os.path.dirname(fasta_path)
    os.chdir(indir)
    cmd = [cmfinder_exe,"-combine", fasta]
    logger.info(f"run cmfinder for {fasta_path} ...")
    p = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=flog)
    flog.close()
            
    os.chdir(cwd)
    if p.returncode == 0:
        logger.info(f"cmfinder run successfully for {fasta_path} .")
        with open(fasta_path + ".checkpoint","w") as f:
            f.write("done")
    else:
        logger.info(f"cmfinder failed for {fasta_path} .")
        with open(fasta_path + ".checkpoint","w") as f:
            f.write("failed")
    return int(p.returncode)

def evaluate_motif():
    flog = open(f"{args.output_directory}/cmsearch.log","wb")
    numbers = {}
    bitscores = {}
    motif_file_name_lookup = {}
    logger.info("evaluate motifs ...")
    for file_name in os.listdir(args.output_directory):
        t = file_name.split(".")[-2]
        if t == "align":
            os.remove(f"{args.output_directory}/{file_name}")
        elif t == "cand":
            os.remove(f"{args.output_directory}/{file_name}")
        elif ".fa.motif." in file_name:
            index = file_name.split(".fa.motif.")[-1]
            if "temp" in index:
                continue 
            motif_file_name_lookup[index] = file_name
            cm = args.output_directory + "/" + index + ".cm"
            cmd = [cmbuild_exe, cm, f"{args.output_directory}/{file_name}"]
            subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=flog)
            fasta = os.path.join(args.output_directory, name + ".90.fa")
            tbl = args.output_directory + "/" + index + ".tbl"
            cmd = [cmsearch_exe,"--toponly", "--tabfile", tbl, cm, fasta]
            subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=flog)
            bitscore = 0
            n = 0
            with open(tbl) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = re.split("\s+",line.strip())
                    score = float(fields[6]) 
                    if score < 10:
                        continue
                    n += 1
                    bitscore += score
            numbers[index] = n
            bitscores[index] = bitscore
    best_index, best_score = -1, -1
    logger.info("bitscores of motif hits:")
    for index in bitscores:
        if best_score < bitscores[index]:
            best_index, best_score = index, bitscores[index]
        print(index, numbers[index], bitscores[index],sep="\t")
    motif = motif_file_name_lookup[best_index]
    logger.info(f"select {motif} as best candidates .")
    logger.info(f"run rscape ...")
    cmd = [rscape_exe, "--cacofold", "--outdir", args.output_directory, motif]
    subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=flog) 
    cm = args.output_directory + "/" + index + ".cm"
    with open(args.output_directory + "/" + "rep.cm","w") as fout:
        fout.write(cm)
    flog.close()

    

def main():
    parser = argparse.ArgumentParser(description='detect conserved structure with cmfinder')
    parser.add_argument('--fasta', '-f', type=str, required=True, help='input RNA sequence in fasta format')
    parser.add_argument('--number', '-n', type=int, default=100, help='number of sequences to sample')
    parser.add_argument('--random-seed', '-rs', type=int, default=666, help='random seed to use')
    parser.add_argument('--output-directory','-od', type=str , help="where to save output, should be absolute path")
    global args
    args = parser.parse_args()
    np.random.seed(args.random_seed)

    if not os.path.exists(args.output_directory):
        logger.info(f"{args.output_directory} does not exists, create it.")
        os.mkdir(args.output_directory)

    global name
    name = args.fasta.split("/")[-1]
    name = name[:name.rfind(".")]
    


    sequences_all, attrs_all = load_fasta(args.fasta)
    logger.info(f"{len(sequences_all)} sequences present in input file .")


    #first clustering sequence to 90% identity
    path90 = os.path.join(args.output_directory,name + ".90.fa")
    if (not os.path.exists(path90)) or rerun_all:
        cd_hits_log = os.path.join(args.output_directory, "cd-hits.log")
        fcd_hits_log = open(cd_hits_log,"w")
        logger.info("clustering input sequences ...")
        cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo/bin/cd-hit-est","-c","0.9", "-r", "0", "-i", args.fasta, "-o", path90,"-d","10000"]    
        subprocess.run(cmd,stdout=fcd_hits_log,stderr=fcd_hits_log)
        fcd_hits_log.close()
    else:
        logger.info(f"already clustered to 90% identity at {path90} .")

    logger.info("pick representative sequence")

    sequences_90, attrs_90 = load_fasta(path90)
    logger.info(f"{len(sequences_90)} sequences after clustering to 90% identity .")
    path90_sampled = os.path.join(args.output_directory, name + ".90.sampled.fa")
    table90 = os.path.join(args.output_directory,name + ".90.txt")
    
    if (not os.path.exists(table90)) or rerun_all:
        clstr90 = os.path.join(args.output_directory,name + ".90.fa.clstr")
        logger.info(f"reformat clustering table as {table90} ...")
        cmd = ["scripts/cd-hit-to-clustering-table-nt.py","-i",clstr90,"-o",table90]
        subprocess.run(cmd)

    clusters = {}
    with open(table90) as f:
        for line in f:
            seq_id, centroid_id = line.strip().split("\t")
            if centroid_id not in clusters:
                clusters[centroid_id] = []
            clusters[centroid_id].append(seq_id)
    centroid_ids = sorted([ (centroid_id, len(seq_ids)) for (centroid_id, seq_ids) in clusters.items()],key=lambda x:-x[-1])

    if (not os.path.exists(path90_sampled)) or rerun_all:
        #if  > n sequence, sampling n sequences
        fout = open(path90_sampled,"w")
        if len(centroid_ids) > args.number:
            logger.info(f"sampling {args.number} in {len(centroid_ids)} sequences for cmfinder search.")
        else:
            logger.info(f"sampling {len(centroid_ids)} sequences for cmfinder search.")
        for i,(centroid_id, size) in enumerate(centroid_ids[:args.number]):
            seq_lengths_in_cluster = []
            seq_ids_in_cluster = []
            for seq_id in clusters[centroid_id]:
                seq_lengths_in_cluster.append(len(sequences_all[seq_id]))
                seq_ids_in_cluster.append(seq_id)
            indices = np.argsort(seq_lengths_in_cluster)
            index = indices[int(len(indices)/2)]
            seq_id = seq_ids_in_cluster[index]
            sequence = sequences_all[seq_id]
            print(f">{i}:{seq_id}", file=fout)
            print(sequence, file=fout)
        fout.close()
        logger.info(f"results saved in {path90_sampled}.")
    else:
        logger.info(f"sampled sequence for cmfinder search presents, see {path90_sampled}.")

    logger.info("run cmfinder ...")
    cmfinder()
    logger.info("run cmsearch ...")
    evaluate_motif() 
    
 
if __name__ == "__main__":
    main()
