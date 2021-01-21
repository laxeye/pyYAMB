#!/usr/bin/env python3
import argparse


def parse_args():
	parser = argparse.ArgumentParser(prog="pyYAMB",
		description=f"pyYAMB metagenome binner")

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible), default 1.")
	general_args.add_argument("-m", "--memory-limit", type=int, default=16,
		help="Memory limit in GB, default 16.")
	input_args = parser.add_argument_group(title="Input files and options")
	input_args.add_argument("-1", "--pe-1", nargs='+',
		help="First (left) paired-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='+',
		help="Second (right) paired-end reads, FASTQ. Space-separated if multiple.")
	input_args.add_argument("-a", "--assembly",
		help="Previously assembled metagenome.")

	asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler",
		default="unicycler", choices=["spades", "megahit"],
		help="Assembler: 'spades' or 'megahit'.")
	asly_args.add_argument("--spades-correction", action="store_true",
		help="Perform short read correction by SPAdes (not recommended).")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")

	args = parser.parse_args()

	'''Check input'''

	return args


def main(args):
	pass

'''Init logger'''


'''Process reads'''

'''Assemble with spades, megahit, ?'''

'''Cut contigs

	fragments = get_fragments(filename, target_length, min_length)

'''

'''Count kmers'''

'''Map with minimap2'''

'''Extract coverage'''

'''Process data'''


if __name__ == '__main__':
	args = parse_args()
	main(args)