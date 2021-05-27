#!/usr/bin/env python3
import logging
from pyyamb.utils import run_external


def map_reads(assembly, reads, target="mapping.sam", threads=1):
	logger = logging.getLogger("main")
	cmd = ["minimap2",
		"-x", "sr", "-t", str(threads),
		"-a", "-o", target, assembly, *reads]

	logger.info("Mapping reads: %s", reads)
	if run_external(cmd).returncode == 0:
		return target
	else:
		logger.error("Unsuccesful mapping.")
		return None
