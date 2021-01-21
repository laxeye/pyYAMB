import subprocess


def run_external(cmd, keep_stdout=False, keep_stderr=False):
	# logger.debug("Running: %s", " ".join(cmd))
	stderr_dest = subprocess.PIPE if keep_stderr else subprocess.DEVNULL
	stdout_dest = subprocess.PIPE if keep_stdout else subprocess.DEVNULL

	r = subprocess.run(cmd, stderr=stderr_dest, stdout=stdout_dest, encoding="utf-8")

	'''if r.returncode != 0:
		logger.error("Non-zero return code of %s", " ".join(cmd))
	'''

	return r


def write_records_to_fasta(records, path):
	with open(path, "w") as f:
		for record in records:
			SeqIO.write(record, f, 'fasta')
		return path
	return None

