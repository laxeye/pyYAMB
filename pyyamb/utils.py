import logging
import subprocess
from Bio import SeqIO


def run_external(cmd, keep_stdout=False, keep_stderr=True):
	logger = logging.getLogger("main")
	logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	try:
		r = subprocess.run(cmd, stderr=stderr_dest, stdout=stdout_dest,
			check=True, encoding="utf-8")
		return r
	except subprocess.CalledProcessError as e:
		logger.error("Error during run of \"%s\"", e.cmd)
		logger.error("stderr message:")
		logger.error(e.stderr)
		raise e


def write_records_to_fasta(records, path):
	'''Return path to FASTA or None'''
	logger = logging.getLogger()
	with open(path, "w") as f:
		for record in records:
			SeqIO.write(record, f, 'fasta')
		return path
		logger.debug("Sequences had been written to: %s", path)
	return None
