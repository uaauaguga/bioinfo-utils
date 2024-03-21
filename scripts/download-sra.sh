#!/bin/bash
dataset=GSE195579-SRP357080
[ -d log/${dataset}/download ] || mkdir -p log/${dataset}/download
[ -d sra/$dataset  ] || mkdir -p sra/$dataset
for run_id in $(cat metadata/${dataset}.txt);do
  if [ -s fastq/$dataset/${run_id}_1.fastq.gz -o -s fastq/$dataset/${run_id}.fastq.gz ];then
    echo "the fastq file of ${run_id} already exists, skip it"
  elif [ -s sra/$dataset/${run_id}.sra ];then 
    echo "the sra file of ${run_id} already exists, skip it"
  else
    echo "downloading ${run_id} ..."
    sem -j 4 "prefetch --resume yes --progress --log-level debug  --max-size 200G --output-file  sra/$dataset/${run_id}.sra ${run_id} > log/${dataset}/download/${run_id}.log 2>&1"
 fi
done
sem --wait
