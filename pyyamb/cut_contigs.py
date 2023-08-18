#!/usr/bin/env python3
import os
import sys
from math import ceil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyyamb.utils import write_records_to_fasta


def get_fragments(filename, target_length=10000, min_length=1000):
	'''Return contig fragments as List[SeqRecord]'''
	fragments = []
	for record in SeqIO.parse(filename, "fasta"):
		seq_length = len(record.seq)
		if seq_length < min_length:
			continue
		if seq_length < 1.5 * target_length:
			fragments.append(record)
		else:
			frag_count = round(seq_length / target_length)
			frag_len = int(ceil(seq_length / frag_count))
			for i in range(frag_count):
				fragments.append(SeqRecord(
					record.seq[frag_len * i: frag_len * (i + 1)],
					id=f"{record.id}_frag_{i}",
					description=""
				))
	return fragments


if __name__ == '__main__':
	filename = sys.argv[1]
	fragments = get_fragments(filename)
	outfilename = f"cut.{os.path.basename(filename)}"
	write_records_to_fasta(fragments, outfilename)
