from pathlib import Path
from sys import argv
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",required=True) #output file
    parser.add_argument("-d","--bam_dir",required=True) #directory containing bamfiles
    parser.add_argument("--min",required=False,default=1.5)
    parser.add_argument("--max",required=False,default=2.5)
    args = parser.parse_args(argv[1:])
    with open(args.file,"w") as outhandle:
        for file in Path(args.bam_dir).glob("*.bam"):
            sample = file.name
            outhandle.write(f"{sample}\t{args.min}\t{args.max}\n")









