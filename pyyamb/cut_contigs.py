#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
from math import ceil
from pyyamb import map


def write_records_to_fasta(records, path):
	with open(path, "w") as f:
		for record in records:
			SeqIO.write(record, f, 'fasta')
		return path
	return None


def get_fragments(filename, target_length=10000, min_length=1000):
	fragments = []
	for record in SeqIO.parse(filename, "fasta"):
		l = len(record.seq)
		if l < min_length:
			continue
		elif l <= target_length:
			fragments.append(record)
		else:
			N = round(l/target_length) # Number of pieces
			L = int(ceil(l/N)) #Length of pieces (except the last one)
			for i in range(N):
				fragments.append(SeqRecord(record.seq[L*i: L*(i+1)],
					id=f"{record.id}_{i}",
					description=""))
	return fragments


if __name__ == '__main__':
	filename = sys.argv[1]
	fragments = get_fragments(filename)
	outfilename = f"cut.{os.path.basename(filename)}"
	write_records_to_fasta(fragments, outfilename)
