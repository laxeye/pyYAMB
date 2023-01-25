#!/usr/bin/env python3
import logging
import os
import os.path
from pyyamb.utils import run_external
import pysam


def map_reads(args):
	logger = logging.getLogger("main")
	if args.single_end:
		# Only first single-end reads file is currently mapped
		reads = [args.single_end[0]]
		prefix, ext = os.path.splitext(os.path.basename(args.single_end[0]))
	if args.pe_1 and args.pe_1:
		# Only one pair of reads is currently mapped
		reads = list(zip(args.pe_1, args.pe_2))[0]
		prefix, ext = os.path.splitext(os.path.basename(args.single_end[0]))

	if ext in ['gz', 'bz2']:
		prefix, _ = os.path.splitext(prefix)
	target = os.path.join(args.output, f"{prefix}.sam")
	cmd = ["minimap2",
		"-x", "sr", "-t", str(args.threads),
		"-a", "-o", target, args.assembly, *reads]

	try:
		logger.info("Mapping reads: %s", ", ".join(reads))
		run_external(cmd)
		return target
	except Exception as e:
		logger.error("Unsuccesful mapping.")
		raise e


def view_mapping_file(args, mapping_sam_file, compress=False):
	logger = logging.getLogger("main")
	logger.info("Converting mapping file")
	prefix = os.path.splitext(mapping_sam_file)
	mapping_bam_file = os.path.join(args.output, f'{prefix}.bam')
	opts = ['-@', str(args.threads), '-F', '0x4', '-o', mapping_bam_file, mapping_sam_file]
	if not compress:
		opts.insert(0, '-u')
	pysam.view(*opts, catch_stdout=False)
	os.remove(mapping_sam_file)

	return mapping_bam_file


def sort_mapping_file(args, mapping_bam_file):
	logger = logging.getLogger("main")
	logger.info("Sorting mapping file")
	prefix = os.path.splitext(mapping_sam_file)
	sorted_mapping_bam_file = os.path.join(args.output, f'{prefix}.sorted.bam')
	'''User-provided memory limit not set here cause it's a memory per thread
	'-m', f'{args.memory_limit}G'
	'''
	pysam.sort('-@', str(args.threads), '-m', '1G', '-o', sorted_mapping_bam_file, mapping_bam_file)
	os.remove(mapping_bam_file)
	logger.info("Indexing mapping file")
	pysam.samtools.index(sorted_mapping_bam_file)

	return sorted_mapping_bam_file
