#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='joint homolog RNA pairs using NNN as separater')
    parser.add_argument('--first-rna','-fr', type=str , required = True, help="query RNA sequence")
    parser.add_argument('--second-rna','-sr', type=str , required = True, help="target RNA sequence")
    parser.add_argument('--output','-o', type=str , required = True, help="output reconstructed probabilities")   
    args = parser.parse_args()
    f1 = open(args.first_rna)  
    f2 = open(args.second_rna)  
    fout = open(args.output,"w")
    for l1,l2 in zip(f1,f2):
        assert l1 == l2
        s1,s2 = next(f1),next(f2)
        print(l1.strip(),file=fout)
        print(s1.strip() + "NNN" + s2.strip(),file=fout)
    f1.close()
    f2.close()
    fout.close()

if __name__ == "__main__":
    main()
