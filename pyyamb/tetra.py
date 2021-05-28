#!/usr/bin/env python3
from Bio import SeqIO
import sys
import regex
import os
import pandas


def compl_DNA(c) -> str:
	'''Return complemetary nucleotide'''
	'''In this script only uppercase nucleotides may appear, otherwise:
	if c.islower():
		c = c.upper()
	'''
	if c == 'A':
		return 'T'
	elif c == 'C':
		return 'G'
	elif c == 'G':
		return 'C'
	elif c == 'T':
		return 'A'
	else:
		return c


def rev_compl_DNA(s) -> str:
	return ''.join(list(map(compl_DNA, s[::-1])))


def make_kmer_list(k=4, nr=True):
	'''Return list of kmers different (if nr) from their reverse complements'''
	acgt = ['A', 'C', 'G', 'T']
	kmers = acgt
	for i in range(k - 1):
		kmers = [f"{x}{y}" for x in kmers for y in acgt]
	if nr is False:
		return kmers
	else:
		nr_kmers = []
		for x in kmers:
			if rev_compl_DNA(x) not in nr_kmers:
				nr_kmers.append(x)
		return nr_kmers


def kmers_freq(records, kmers):
	'''Yield kmer frequence in FASTA file'''
	patterns = [(i, regex.compile(i)) for i in kmers]
	klen = len(kmers[0])
	for record in records:
		d = {}
		seq_len = len(record.seq)
		seq = str(record.seq)
		for (i, j) in patterns:
			d[i] = len(regex.findall(j, seq, overlapped=True))
		seq_rc = str(record.seq.reverse_complement())
		for (i, j) in patterns:
			d[i] = d.get(i, 0) + len(regex.findall(j, seq_rc, overlapped=True))
		yield [record.id, seq_len] + [1000 * d.get(i, 0) / (2 * (seq_len - klen)) for i in kmers]


def kmer_freq_table(filename, k_len=4):
	kmer_list = make_kmer_list(k_len, nr=True)
	return pandas.DataFrame.from_records(
		kmers_freq(SeqIO.parse(filename, "fasta"), kmer_list),
		columns=['fragment', 'length'] + kmer_list,
		index='fragment')


if __name__ == '__main__':
	outfilename = f"tn.{os.path.basename(sys.argv[1])}.csv"
	kmer_freq_table(sys.argv[1]).to_csv(outfilename)
