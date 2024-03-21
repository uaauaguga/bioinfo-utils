#!/bin/bash
dataset=GSE195579-SRP357080
mkdir -p log/${dataset}/dump fastq/$dataset
for run_id in $(cat metadata/$dataset.txt );do
if [ -s fastq/$dataset/${run_id}_1.fastq.gz -o -s fastq/$dataset/${run_id}.fastq.gz ];then
    echo "the fastq file of ${run_id} already exists, skip it"
  elif [ ! -s sra/$dataset/${run_id}.sra ];then 
    echo "the sra file of ${run_id} do not exists, skip it"
  else
    echo "dump ${run_id} ..."
    #sem -j 10 "fastq-dump --gzip --split-3 sra/${dataset}/${run_id}.sra -O fastq/$dataset  > log/${dataset}/dump/${run_id}.txt 2>&1"
 fi
done
#sem --wait
