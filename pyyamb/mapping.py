#!/usr/bin/env python3

import argparse
import typing
import subprocess
import os


def map_reads(args, assembly, reads):
	mapping = os.path.join(args.output_dir, "mapping.sam")
	if args.mapper == 'minimap2':
		result = map_minimap(args, assembly, reads, mapping)
	'''
	elif args.mapper == 'bwa':
		map_bwa()
	elif args.mapper =='bowtie2':
		map_bowtie() 
	'''
	if result.returncode == 0:
		return mapping
	else:
		return None


def map_minimap(args, assembly, reads, mapping):
	threads = args.threads if args.threads else 1
	cmd = ['minimap2', '-t', str(args.threads), '-a', '-x', 'sr', '-o', mapping]
	if isinstance(reads, (list, tuple, set)):
		cmd += list(reads)
	else:
		cmd.append(reads)
	try:
		return subprocess.run(cmd)
	except Exception as e:
		raise e
	else:
		return None


def parse_args():
	parser = argparse.ArgumentParser(description="pyYAMB mapping module")

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--mapper", help="Mapping software", choices=["minimap2", "bwa", "bowtie2"], default="minimap2")
	general_args.add_argument("-o", "--output-dir", help="Output directory")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible)")
	input_args = parser.add_argument_group(title="Input options", description="")
	input_args.add_argument("-a", "--assembly", help="Metagenomic assembly")
	input_args.add_argument("-1", "--pe-1", help="FASTQ file with first (left) paired-end reads")
	input_args.add_argument("-2", "--pe-2", help="FASTQ file with second (right) paired-end reads")
	input_args.add_argument("-s", "--single-end", help="FASTQ file with second (right) paired-end reads")

	args = parser.parse_args()

	if not args.output_dir:
		args.output_dir = os.getcwd()

	return args


def main():
	args = parse_args()
	reads = []
	if os.path.isfile(args.pe_1) and os.path.isfile(args.pe_2):
		reads += (os.path.abspath(args.pe_1), os.path.abspath(args.pe_2))
	if os.path.isfile(args.single_end):
		reads += os.path.abspath(args.single_end)
	print(map_reads(args, args.assembly, reads))
	if len(reads) == 0:
		raise FileNotFoundError("Input files not found!")


if __name__ == '__main__':
	main()
