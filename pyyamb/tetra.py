#!/usr/bin/env python3
from Bio import SeqIO
import sys
import regex
import os


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


def rev_compl_DNA(s) -> str :
	return ''.join(list(map(compl_DNA, s[::-1])))


def make_kmer_list(k=4, nr=True):
	'''Return list of kmers different (if nr) from their reverse complements'''
	acgt = ['A','C','G','T']
	kmers = acgt
	for i in range(k-1):
		kmers = [f"{x}{y}" for x in kmers for y in acgt]
	if nr is False:
		return kmers
	else:
		nr_kmers = []
		for x in kmers:
			if rev_compl_DNA(x) not in nr_kmers:
				nr_kmers.append(x)
		return nr_kmers

	
def kmers_freq(seqlist, kmers):
	'''Yield kmer frequence in DNA sequence list'''
	patterns = [(i, regex.compile(i)) for i in kmers]
	klen = len(kmers[0])
	for (id, seq) in seqlist:
		d={}
		l = len(seq)
		for (i,j) in patterns:
			d[i] = len(regex.findall(j, str(seq), overlapped=True))
		for (i,j) in patterns:
			d[i] = d.get(i, 0) + len(regex.findall(j, str(seq.reverse_complement()), overlapped=True))
		p = "\t".join([str(1000*d.get(i,0) / (2*(l - klen))) for i in tnlist])
		yield f"{id}\t{p}"


def main(filename):
	outfilename = f"tn.{os.path.basename(filename)}.tsv"
	g = ((record.id, record.seq) for record in SeqIO.parse(filename, "fasta"))
	kmer_list = make_kmer_list(k=4, nr=True)
	with open(outfilename, 'w') as OF:
		for lane in kmers_freq(g, kmer_list):
			print(lane, file=OF)


if __name__ == '__main__':
	main(sys.argv[1])

