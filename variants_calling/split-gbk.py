#!/usr/bin/env python

'''
This is as script to split genbank files into separated sequences
Here, we only keep the chromosome sequence
Later we will merge in one script and use an option instead
'''

import sys
from Bio import SeqIO

def split_gbk(file):
    with open(file,'r') as input_handle:
        seqs=SeqIO.parse(input_handle,'genbank')
        for record in seqs:
            header=record.id
            with open(record.id+'.gbk','w') as output_handle:
                SeqIO.write(record, output_handle,'genbank')

def main():
    if len(sys.argv) == 2:
        file = sys.argv[1]
        split_gbk(file)
    else:
        print("ERROR: split-gbk.py takes one argument.")

if __name__ == "__main__":
    main()

        
        
