#!/usr/bin/env python
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("concatenate table")


def main():
    parser = argparse.ArgumentParser(description='concatenate table in a directory')
    parser.add_argument('--input-directory', '-id', required=True, help="input intervals")
    parser.add_argument('--output','-o', required=True, help="output intervals")
    parser.add_argument('--prefix','-p', action = "store_true", help="whether appened a prefix")
    parser.add_argument('--input-suffix','-is', help="only consider input with this suffix")
    args = parser.parse_args()

    logger.info(f"saving results to {args.output} ...")
    fout = open(args.output,"w")
    for txt in tqdm(sorted(os.listdir(args.input_directory))):
        if args.input_suffix and not txt.endswith(args.input_suffix):
            continue
        sample_id = txt[:-4]
        path = os.path.join(args.input_directory,txt)
        try:
            with open(path) as f:
                for line in f:
                    if args.prefix and (not line.startswith("#")):
                        line = sample_id + ":" + line
                    fout.write(line)
        except:
            logger.info(f"some thing goes wrong with {sample_id}, skip it .")
    fout.close()
    logger.info("all done .")


if __name__ == "__main__":
    main()
    
