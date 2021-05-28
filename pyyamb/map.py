#!/usr/bin/env python3
import logging
import os
from pyyamb.utils import run_external


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
