#!/usr/bin/env python3
import logging
import os
from pyyamb.utils import run_external
import pysam


def map_reads(args):
	logger = logging.getLogger("main")
	target = os.path.join(args.output, "mapping.sam")
	if args.single_end:
		# Only first single-end reads file is currently mapped
		reads = [args.single_end[0]]
	if args.pe_1 and args.pe_1:
		# Only one pair of reads is currently mapped
		reads = list(zip(args.pe_1, args.pe_2))[0]

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
	mapping_bam_file = os.path.join(args.output, 'mapping.bam')
	opts = ['-@', str(args.threads), '-F', '0x4', '-o', mapping_bam_file, mapping_sam_file]
	if compress:
		opts.insert(0, '-u')
	pysam.view(*opts, catch_stdout=False)
	os.remove(mapping_sam_file)

	return mapping_bam_file


def sort_mapping_file(args, mapping_bam_file):
	logger = logging.getLogger("main")
	logger.info("Sorting mapping file")
	sorted_mapping_bam_file = os.path.join(args.output, 'mapping.sorted.bam')
	pysam.sort('-@', str(args.threads), '-o', sorted_mapping_bam_file, mapping_bam_file)
	os.remove(mapping_bam_file)
	logger.info("Indexing mapping file")
	pysam.samtools.index(sorted_mapping_bam_file)

	return sorted_mapping_bam_file
