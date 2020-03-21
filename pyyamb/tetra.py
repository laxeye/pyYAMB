#!/usr/bin/env python3
from Bio import SeqIO
import sys
import regex

#Upper case only :(
def compl(c):
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

def rev_compl(s) -> str :
	return ''.join(list(map(compl, s[::-1])))

#Create a list of tetranucleotides different from their revcomp
def make_tn_list():
	nr_tn = []
	acgt = ['A','C','G','T']
	tnlist = [f"{a}{b}{c}{d}" for a in acgt for b in acgt for c in acgt for d in acgt]
	for x in tnlist:
		if rev_compl(x) not in nr_tn:
			nr_tn.append(x)
	return nr_tn
	
def tn_freq(seqlist):
	tnlist = make_tn_list()
	patterns = [(i, regex.compile(i)) for i in tnlist]
	
	for (id, seq) in seqlist:
		d={}
		l = len(seq)
		for (i,j) in patterns:
			d[i] = len(regex.findall(j, str(seq), overlapped=True))
		for (i,j) in patterns:
			d[i] = d.get(i, 0) + len(regex.findall(j, str(seq.reverse_complement()), overlapped=True))
		p = "\t".join([str(1000*d.get(i,0)/(2*l - 6)) for i in tnlist])
		yield f"{id}\t{p}"

def main(infilename):
	outfilename = f"{infilename}.csv"
	g = ((record.id, record.seq) for record in SeqIO.parse(infilename, "fasta"))
	with open(outfilename, 'w') as OF:
		for lane in tn_freq(g):
			print(lane, file=OF)

if __name__ == '__main__':
	main(sys.argv[1])

