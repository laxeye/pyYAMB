#!/usr/bin/env python3
import argparse
from pyyamb import cutcontigs
from pyyamb import mapping
import typing
import os


def parse_args():
	parser = argparse.ArgumentParser(prog="pyYAMB",
		description=f"pyYAMB metagenome binner")

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", help="Task of pipeline",
		choices=["assemble", "cut", "map", "bin", "checkm"], required=True)
	general_args.add_argument("--mapper", help="Mapping software", choices=["minimap2", "bwa", "bowtie2"], default="minimap2")
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

	parser.add_argument("--minlength", help="minimum contig length", default=1000, type=int)
	parser.add_argument("--fraglength", help="Target length of contig fragments", default=10000, type=int)

	args = parser.parse_args()

	'''Check input'''

	return args


def main(args):
	'''Init logger'''

	'''Process reads'''

	'''Assemble with spades, megahit, ?'''

	if args.task == "cut":
		try:
			args.assembly = os.path.abspath(args.assembly)
		except Exception as e:
			raise e
		cutcontigs.write_cut_contigs(args.assembly, args.fraglength, args.minlength)

	if args.task == "map":
		if os.path.isfile(args.pe_1) and os.path.isfile(args.pe_2):
			reads = (os.path.abspath(args.pe_1), os.path.abspath(args.pe_2))
			if mapping.map_reads(args, args.assembly, reads):
				print("Reads mapped!")
		else:
			raise FileNotFoundError(f"Input files {args.pe_1} and {args.pe_2} not found!")



	'''Cut contigs

		fragments = get_fragments(filename, target_length, min_length)

	'''

	'''Count kmers'''

	'''Map with minimap2'''

	'''Extract coverage'''

	'''Process data'''


if __name__ == '__main__':
	args = parse_args()
	main()
