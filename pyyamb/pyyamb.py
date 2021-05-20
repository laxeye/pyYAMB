#!/usr/bin/env python3

import argparse
from pyyamb import cutcontigs
from pyyamb import mapping
import typing
import os


def parse_args():
	parser = argparse.ArgumentParser(description="pyYAMB")

	# General options
	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", help="Task of pipeline",
		choices=["assemble", "cut", "map", "bin", "checkm"], required=True)
	general_args.add_argument("--mapper", help="Mapping software", choices=["minimap2", "bwa", "bowtie2"], default="minimap2")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible)")

	input_args = parser.add_argument_group(title="Input options", description="")
	input_args.add_argument("-a", "--assembly", help="Metagenomic assembly")
	input_args.add_argument("-1", "--pe-1", help="FASTQ file with first (left) paired-end reads")
	input_args.add_argument("-2", "--pe-2", help="FASTQ file with second (right) paired-end reads")

	parser.add_argument("--minlength", help="minimum contig length", default=1000, type=int)
	parser.add_argument("--fraglength", help="Target length of contig fragments", default=10000, type=int)

	return parser.parse_args()


def main():
	args = parse_args()
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


if __name__ == '__main__':
	main()
