#!/usr/bin/env python
import argparse
from collections import defaultdict
from multiprocessing import Pool
import os
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('build msa')


cd_hit_est = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/cd-hit-est"
mafft = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/mafft"
hmmbuild = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/hmmbuild"
hmmsearch = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/hmmsearch"
nhmmer = "/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/nhmmer"

def build_msa(cluster_id):
    output_directory = os.path.join(args.output_directory,cluster_id[-3:],cluster_id)
    fasta = os.path.join(output_directory, cluster_id + ".fa")
    ff = open(fasta,"w")
    for header, sequence in groupped_sequences[cluster_id]:
        ff.write(header)
        ff.write(sequence)
    ff.close()
    # clustering to 80% identity
    logger.info(f"clustering {cluster_id} to 80% identity ...")
    fasta80 = os.path.join(output_directory, cluster_id + ".80.fa")
    cmd = [cd_hit_est, "-i", fasta, "-o", fasta80, "-r", "0", "-c", "0.8"]
    fcd_hit_log = open(os.path.join(output_directory, cluster_id + ".cd-hit.log"),"w")
    subprocess.run(cmd,stdout=fcd_hit_log,stderr=fcd_hit_log)
    fcd_hit_log.close()

    # sampling up to 100 sequences
    fasta80_sampled = os.path.join(output_directory, cluster_id + ".80.sampled.fa")
    sequences = {}
    with open(fasta80) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip()
    with open(fasta80_sampled,"w") as ff:
        for i,seq_id in enumerate(sequences):
            if i >= 100:
                break
            print(f">{seq_id}",file=ff)
            print(sequences[seq_id] + "\n",file=ff)
            
    sampled_msa = os.path.join(output_directory, cluster_id + ".80.sampled.aligned.fa")    
    fmsa = open(sampled_msa,"w")
    # run mafft
    logger.info(f"run mafft on {min(100,len(sequences))} {cluster_id} sequences ...")
    cmd = [mafft,"--maxiterate","1000","--localpair",fasta80_sampled]
    fmafft_log = open(os.path.join(output_directory, cluster_id + ".mafft.log"),"w")
    subprocess.run(cmd,stdout=fmsa,stderr=fmafft_log)
    fmsa.close()
    fmafft_log.close()

    # build hmmprofile
    fhmmbuild_log = open(os.path.join(output_directory, cluster_id + ".sampled.hmmbuild.log"),"w")
    sampled_hmm = os.path.join(output_directory, cluster_id + ".80.sampled.hmm")    
    cmd = [hmmbuild, "-n", cluster_id + "-sampled", sampled_hmm, sampled_msa]
    subprocess.run(cmd,stdout=fhmmbuild_log,stderr=fhmmbuild_log)
    fhmmbuild_log.close()

    # build full msa
    logger.info(f"build full msa for {cluster_id} ...")
    full_msa = os.path.join(output_directory, cluster_id + ".aligned.stk")
    tbl = os.path.join(output_directory, cluster_id + ".tbl")
    #cmd = [hmmsearch,"--tblout",tbl,"-A",full_msa,sampled_hmm, fasta]
    cmd = [nhmmer,"--watson","--tblout",tbl,"-A",full_msa,sampled_hmm, fasta]
    subprocess.run(cmd,stdout=subprocess.DEVNULL)

    # build full hmm profile
    fhmmbuild_log = open(os.path.join(output_directory, cluster_id + ".full.hmmbuild.log"),"w")
    full_hmm = os.path.join(output_directory, cluster_id + ".hmm")    
    cmd = [hmmbuild, "-n", cluster_id + "-full", full_hmm, full_msa]
    subprocess.run(cmd, stdout=fhmmbuild_log,stderr=fhmmbuild_log)
    fhmmbuild_log.close()
    logger.info(f"hmm profile saved in {full_hmm} .")
    return 0

def main():
    global args
    parser = argparse.ArgumentParser(description='build msa of groupped sequences')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences in fasta format')
    parser.add_argument('--output-directory','-od',type=str, required=True, help="output directory")
    parser.add_argument('--jobs','-j',type=int, default=8, help='number of jobs')
    args = parser.parse_args()
    global groupped_sequences
    groupped_sequences = defaultdict(list)
    logger.info("load groupped sequences ...")
    with open(args.input) as f:
        for header in f:
            sequence = next(f)
            cluster_id = header.split(" ")[2]
            groupped_sequences[cluster_id].append((header,sequence))

    if not os.path.exists(args.output_directory):
        try:
            os.mkdir(args.output_directory)
        except:
            pass
    pool = Pool(args.jobs) 
    workers = []
    for cluster_id in groupped_sequences:
        if len(groupped_sequences[cluster_id]) < 10:
            continue
        odp = os.path.join(args.output_directory,cluster_id[-3:])
        if not os.path.exists(odp):
            os.mkdir(odp)
        od = os.path.join(args.output_directory,cluster_id[-3:],cluster_id)
        if not os.path.exists(od):
            try:
                os.mkdir(od)
            except:
                continue
        workers.append(pool.apply_async(func=build_msa, args=(cluster_id,))) 
    for worker in workers:
        code = worker.get() 



if __name__ == "__main__":
    main()
