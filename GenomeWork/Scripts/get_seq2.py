#!/usr/bin/python
from Bio import SeqIO
import sys                                                                          

"""
usage: ./get_seq2.py input.fasta wanted_taxa.txt output.fasta
"""

# give a fasta file with lots of sequences, a taxa file with taxa on different lines that you want seq from
# and an output fasta file name


fasta_file = SeqIO.parse(open(sys.argv[1], 'r'), "fasta")	# Input fasta file
wanted_file = open(sys.argv[2], 'r')						# Input interesting sequence IDs, one per line
result_file = open(sys.argv[3], 'w')						# Output fasta file

wanted = []													# make empty list and fill with taxa from your wanted file
with wanted_file as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.append(line)
            print(line)


with result_file as f:										# loop through list and if the name in the taxa list is in the sequence
    for seq in fasta_file:									# description it will put it into the new output file
    	for taxa in wanted:
    		if taxa in seq.description:
	            SeqIO.write(seq, f, "fasta")
