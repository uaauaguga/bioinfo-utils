#!/usr/bin/env python
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='split text files by group name, make sure adjacent lines belones to the same group')
    parser.add_argument('--input', '-i', type=str, required=True, help='input file in bed format')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='directory to save bed files')
    parser.add_argument('--extractor','-e', type=str, default="lambda x:x.split(':')[0]", help="function to extract group name from each line")
    parser.add_argument('--renamer','-r', type=str, default="lambda x,n:x", help="how to name the files")
    parser.add_argument('--with-header','-wh', action="store_true", help="whether input data contains header")
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    i = 0
    n = 0
    suffix = args.input.split(".")[-1] 
    extractor = eval(args.extractor)
    renamer = eval(args.renamer)
    last_group_id = ""
    fout = None
    header = None
    n = 0
    with open(args.input) as f:
        if args.with_header:
            header = next(f)
        for line in f:
            group_id = extractor(line)
            name = renamer(group_id,n)
            if group_id != last_group_id:
                if fout is not None:
                    fout.close() 
                fout = open(os.path.join(args.output_directory, name + "." + suffix),"w")
                if header is not None:
                    fout.write(header)
                n += 1
            fout.write(line)            
            last_group_id = group_id
    fout.close()

if __name__ == "__main__":
    main()
