#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
from math import ceil

def get_cut_contigs(filename, target_length, min_length):
	result = []
	f = open(f"cut.{os.path.basename(filename)}", "w")
	for record in SeqIO.parse(filename, "fasta"):
		l = len(record.seq)
		if l < min_length:
			continue
		elif l <= target_length:
			result.append(record)
		else:
			N = round(l/target_length) # Number of pieces
			L = int(ceil(l/N)) #Length of pieces (except the last one)
			i = 0
			while i < N:
				result.append(SeqRecord(record.seq[L*i: L*(i+1)],id=f"{record.id}_{i}"))
				i += 1
	return result


def write_cut_contigs(filename, target_length, min_length):
	'''
	Ignore short contigs and cut long contigs to same-length pieces.
	'''
	f = open(f"cut.{os.path.basename(filename)}", "w")
	for record in SeqIO.parse(filename, "fasta"):
		l = len(record.seq)
		if l < min_length:
			continue
		elif l <= target_length:
			SeqIO.write(record, f, "fasta")
		else:
			N = round(l/target_length) # Number of pieces
			L = int(ceil(l/N)) #Length of pieces (except the last one)
			i = 0
			while i < N:
				SeqIO.write(SeqRecord(record.seq[L*i: L*(i+1)],id=f"{record.id}_[{i}]", description=''), f, 'fasta')
				i += 1
	f.close()

def main():
	write_cut_contigs(sys.argv[1], 10000, 1000)

if __name__ == '__main__':
	main()
