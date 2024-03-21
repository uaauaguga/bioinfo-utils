#!/usr/bin/env python
import argparse 
import logging
from itertools import product
import numpy as np
import json
from sklearn.metrics import roc_curve,auc
from lightgbm import LGBMClassifier
import pickle
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("evaluate predictive power of k-mer profile")


def load_fasta(path):
    sequences = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                seq_id  = line.strip()[1:].split(" ")[0]
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.strip().upper().replace("T","U")
    seq_ids = list(sequences.keys())
    for seq_id in seq_ids:
        if "N" in sequences[seq_id]:
            del sequences[seq_id]
    return sequences


def get_kmer_profile(sequences):
    seq_ids = []
    kmer_profiles = []
    token_size = len(tokens_lut)
    omited_tokens = set()
    for seq_id in sequences:
        seq_ids.append(seq_id)
        x = np.zeros(token_size,dtype=int)
        sequence = sequences[seq_id]
        for i in range(len(sequence)-k):
            token = sequence[i:i+k]
            if token not in tokens_lut:
                omited_tokens.add(token)
                continue
            x[tokens_lut[token]] += 1
        kmer_profiles.append(100*x/len(sequence))
    kmer_profiles = np.array(kmer_profiles)
    logger.warning("These tokens were omited:")
    logger.warning(",".join(omited_tokens))
    return seq_ids, kmer_profiles


def main():
    parser = argparse.ArgumentParser(description='sampling negative sequences')
    parser.add_argument('--train-positive', type=str, required=True, help='positive instance for training')
    parser.add_argument('--train-negative', type=str, required=True, help="negative instance for training")
    parser.add_argument('--test-positive', type=str, required=True, help='positive instance for testing')
    parser.add_argument('--test-negative', type=str, required=True, help="negative instance for testing")
    parser.add_argument('--k-mer', '-k', type=int, default=3, help="k-mer size to use")
    parser.add_argument('--seed','-s',type=int,default=666,help="Random seed for sampling")
    parser.add_argument('--n-jobs','-j',type=int,default=4,help="Number of worker to use for the classifier")
    parser.add_argument('--model','-m',type=str,help="Where to save the model in pickle format")

    args = parser.parse_args()
    tokens = []
    for token in product("ACGU",repeat=args.k_mer):
        tokens.append("".join(token))
    global k
    global tokens_lut
    tokens_lut = {}
    k = args.k_mer
    for i, token in enumerate(sorted(tokens)):
        tokens_lut[token] = i
    logger.info("Load positive training sequences ...")
    train_positive_sequences = load_fasta(args.train_positive)
    logger.info(f"Calculate {k}-mer profile ...")
    train_positive_seq_ids, train_positive_kmer_profiles = get_kmer_profile(train_positive_sequences)
   
    
    logger.info("Load negative training sequences ...")
    train_negative_sequences = load_fasta(args.train_negative)
    logger.info(f"Calculate {k}-mer profile ...")
    train_negative_seq_ids, train_negative_kmer_profiles = get_kmer_profile(train_negative_sequences)
   
    y_train = [1]*len(train_positive_seq_ids) + [0]*len(train_negative_seq_ids)
    y_train = np.array(y_train) 
    X_train = np.concatenate([train_positive_kmer_profiles,train_negative_kmer_profiles],axis=0)
    
    train_idx = np.arange(y_train.shape[0]) 
    np.random.shuffle(train_idx)
    y_train = y_train[train_idx]
    X_train = X_train[train_idx,:]

    logger.info("Fitting GBDT classifier ...")
    clf = LGBMClassifier(boosting_type='gbdt',n_estimators=1000,random_state=args.seed, n_jobs= args.n_jobs,verbose=1)
    clf.fit(X_train,y_train)
  
 
    logger.info("Load positive testing sequences ...")
    test_positive_sequences = load_fasta(args.test_positive)
    logger.info(f"Calculate {k}-mer profile ...")
    test_positive_seq_ids, test_positive_kmer_profiles = get_kmer_profile(test_positive_sequences)
  
    logger.info("Load negative testing sequences ...")
    test_negative_sequences = load_fasta(args.test_negative)
    logger.info(f"Calculate {k}-mer profile ...")
    test_negative_seq_ids, test_negative_kmer_profiles = get_kmer_profile(test_negative_sequences)

    y_test = [1]*len(test_positive_seq_ids) + [0]*len(test_negative_seq_ids)
    y_test = np.array(y_test) 
    X_test = np.concatenate([test_positive_kmer_profiles,test_negative_kmer_profiles],axis=0)
    
    test_idx = np.arange(y_test.shape[0]) 
    np.random.shuffle(test_idx)
    y_test = y_test[test_idx]
    X_test = X_test[test_idx,:]

    logger.info("Make prediction on testing instances ...")
    y_pred_proba = clf.predict_proba(X_test)[:,1]
    fpr, tpr,_ = roc_curve(y_test, y_pred_proba)
    AUROC = auc(fpr, tpr)
    logger.info(f"The testing AUROC score is: {AUROC} ")

    logger.info(f"Saving the model to {args.model} ...")
    with open(args.model,"wb") as f:
        pickle.dump(clf,f)
    logger.info("All done .")

if __name__ == "__main__":
    main()
