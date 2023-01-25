#!/usr/bin/env python3
import argparse
import logging
import os
import pandas
import seaborn
import hdbscan
import pysam
from Bio import SeqIO
from matplotlib import pyplot
from numpy import mean
from sklearn.manifold import TSNE
from pyyamb.cut_contigs import get_fragments
from pyyamb.map import map_reads
from pyyamb.map import view_mapping_file
from pyyamb.map import sort_mapping_file
from pyyamb.utils import write_records_to_fasta
from pyyamb.tetra import kmer_freq_table


def parse_args():
	logger = logging.getLogger()
	parser = argparse.ArgumentParser(prog="pyYAMB",
		description="pyYAMB metagenome binner",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", required=True,
		choices=["all", "cut", "tetra", "clustering", "map", "make_bins"],
		help="Task of pipeline: cut (discard short contigs and cut longer), "
		+ "map (map reads or process mapping file)")
	general_args.add_argument("-o", "--output", type=str, required=True,
		help="Output directory")
	'''general_args.add_argument("--mapper", help="Mapping software, currently minimap2 only",
		choices=["minimap2"], default="minimap2")'''
	general_args.add_argument("--min-length", default=1000, type=int,
		help="Minimum contig length")
	general_args.add_argument("--fragment-length", default=10000, type=int,
		help="Target length of contig fragments")
	general_args.add_argument("--perplexity", default=50, type=int,
		help="Perplexity parameter for tSNE.")
	general_args.add_argument("--k-len", default=4, type=int,
		help="Length of k-mer to calculate their frequencies.")
	general_args.add_argument("--log-norm-coverage", action='store_true',
		help="Perform log-normalization of coverage.")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible).")
	general_args.add_argument("-m", "--memory-limit", type=int, default=16,
		help="Memory limit in GB.")

	input_args = parser.add_argument_group(title="Input files and options")
	input_args.add_argument("-1", "--pe-1", nargs='*',
		help="First (left) paired-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='*',
		help="Second (right) paired-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-s", "--single-end", nargs='*',
		help="Sinle-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-i", "--assembly",
		help="Previously assembled metagenome.")
	input_args.add_argument("--kmers-data",
		help="Previously calculated kmer-freqs.")
	input_args.add_argument("--mapping-file", nargs='*',
		help="Sorted and indexed BAM-file(s). Space-separated if multiple.")
	input_args.add_argument("--coverage-data",
		help="Coverage depth of fragments in comma-separated file.")

	'''asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler",
		default="spades", choices=["spades", "megahit"],
		help="Assembler: 'spades' or 'megahit'.")
	asly_args.add_argument("--spades-correction", action="store_true",
		help="Perform short read correction by SPAdes (not recommended).")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")'''

	args = parser.parse_args()

	args.output = os.path.abspath(args.output)

	if not os.path.isdir(args.output):
		try:
			os.makedirs(args.output)
		except Exception as e:
			logger.error("Failed to create %s", args.output)
			raise e

	'''Check input'''
	for x in (args.assembly, args.kmers_data):
		if x:
			x = os.path.abspath(x)
			if not os.path.isfile(x):
				logger.error("File not found: %s", x)
				raise FileNotFoundError(x)

	if args.single_end:
		args.single_end = list(map(os.path.abspath, args.single_end))
		for x in args.single_end:
			if not os.path.isfile(x):
				logger.error("File not found: %s", x)
				raise FileNotFoundError(x)

	if args.pe_1:
		args.pe_1 = list(map(os.path.abspath, args.pe_1))
		args.pe_2 = list(map(os.path.abspath, args.pe_2))
		for f in args.pe_1 + args.pe_2:
			if not os.path.isfile(f):
				raise FileNotFoundError(f"Input file {f} not found!")

	if args.mapping_file:
		args.mapping_file = list(map(os.path.abspath, args.mapping_file))
		for m_file in args.mapping_file:
			if not os.path.isfile(m_file):
				logger.error("File not found: %s", m_file)
				raise FileNotFoundError(f"{m_file} not found!")

	return args


def create_logger():
	'''Create logger'''
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	return logger


def extract_coverage(bamfile, fasta):
	'''Extract coverage from sorted and indexed BAM-file for every contig in FASTA-file

	Args:
		bamfile (path): Path to sorted and indexed mapping file in BAM format
		fasta (path): Path to metagenome fragments in FASTA format

	Returns:
		pandas.DataFrame
	'''
	sample_name, _ = os.path.splitext(os.path.basename(bamfile))
	contig_list = [record.id for record in SeqIO.parse(fasta, 'fasta')]
	with pysam.AlignmentFile(bamfile) as bam_handle:
		covs = [mean(bam_handle.count_coverage(contig=contig)) for contig in contig_list]
	return pandas.DataFrame.from_records(zip(contig_list, covs), index='fragment', columns=['fragment', sample_name])


def make_nice_bam(args):
	sam_file = map_reads(args)
	bam_file = view_mapping_file(args, sam_file, compress=False)
	return sort_mapping_file(args, bam_file)


def main():
	logger = create_logger()
	args = parse_args()
	logger.info("Analysis started")
	mg_data = pandas.DataFrame()

	'''Process reads'''

	'''Assemble with spades or megahit'''

	'''Cut contigs'''
	if args.task == "cut" or args.task == "all":
		try:
			fragments = get_fragments(args.assembly, args.fragment_length, args.min_length)
			args.assembly = write_records_to_fasta(fragments, os.path.join(args.output, 'fragments.fasta'))
			logger.info("Contigs fragmented.")
		except Exception as e:
			logger.error("Error during contig fragmentation.")
			raise e

	'''Count kmers'''
	if args.task == "tetra" or args.task == "all":
		try:
			logger.info("Counting kmers")
			mg_data = kmer_freq_table(args.assembly, args.k_len, args.threads)
			kmers_data = os.path.join(args.output, "kmers.csv")
			mg_data.to_csv(kmers_data)
			logger.info('Wrote k-mer frequencies to %s.', kmers_data)
		except Exception as e:
			logger.error("Error during kmer frequency calculation.")
			raise e

	'''Map reads with minimap2'''
	if args.task == "map" or (args.task == 'all' and (args.pe_1 or args.single_end)):
		sorted_bam_files = args.mapping_file if args.mapping_file else [make_nice_bam(args)]
		logger.info("Extracting coverage")
		cov_data = extract_coverage(sorted_bam_files[0], args.assembly)
		if len(sorted_bam_files) > 1:
			for m_file in sorted_bam_files[1:]:
				cov_data = cov_data.merge(extract_coverage(m_file, args.assembly), how='outer')
		cov_data.to_csv(os.path.join(args.output, "coverage.csv"))
		logger.info("Processing of mapping file finished")

	'''Cluster data with HDBSCAN'''
	if args.task == "clustering":
		'''Read k-mers from disk. Frequencies, not z-scores'''
		if args.kmers_data and os.path.isfile(args.kmers_data):
			mg_data = pandas.read_csv(args.kmers_data, index_col=0)
			logger.info('Read k-mer frequencies from %s.', args.kmers_data)
		else:
			logger.error("File with k-mers not found")
			raise FileNotFoundError(args.kmers_data)

		'''Read coverage depth of fragments from disk. Not z-scores'''
		if args.coverage_data and os.path.isfile(args.coverage_data):
			cov_data = pandas.read_csv(args.coverage_data, index_col=0)
			logger.info('Read fragment coverage depth from %s.', args.coverage_data)
		else:
			logger.error("File fragment coverage depth not found")
			raise FileNotFoundError(args.coverage_data)

	if args.task == "clustering" or args.task == "all":
		'''Merge data and make Z-score normalization'''
		mg_data = mg_data.merge(cov_data, how='outer', on='fragment')
		mg_data.iloc[:, 2:] = mg_data.iloc[:, 2:].apply(
			lambda x: (x - x.mean()) / x.std(), axis=0)

		logger.info('tSNE data reduction.')
		tsne_pca = TSNE(init='pca', perplexity=args.perplexity).fit_transform(mg_data.iloc[:, 2:])
		dfTSNE = pandas.DataFrame.join(
			pandas.DataFrame(tsne_pca, columns=['tsne1', 'tsne2'], index=mg_data.index),
			mg_data['length']
		)

		logger.info('HDBSCAN data clustering.')
		clusterer = hdbscan.HDBSCAN(min_cluster_size=25)
		cluster_labels = clusterer.fit_predict(dfTSNE[['tsne1', 'tsne2']])
		dfTSNE = dfTSNE.join(pandas.DataFrame(
			cluster_labels, index=mg_data.index, columns=['cluster']
		))
		dfTSNE.to_csv(os.path.join(args.output, f"pyyamb.perplexity_{args.perplexity}.csv"))
		clusters = set(dfTSNE['cluster'])
		logger.info("HDBSCAN found %s clusters", len(clusters))

		logger.info("Producing graphics")
		pal = seaborn.color_palette(palette="Set3", n_colors=len(clusters))
		seaborn.relplot(
			data=dfTSNE, x='tsne1', y='tsne2', size='length',
			alpha=0.2, hue=dfTSNE["cluster"].astype("category"),
			palette=pal
		)
		pyplot.savefig(os.path.join(args.output, f"pyyamb.perplexity_{args.perplexity}.png"), dpi=300)
		pyplot.savefig(os.path.join(args.output, f"pyyamb.perplexity_{args.perplexity}.svg"))

	if args.task in ["make_bins", "clustering", "all"]:
		dfTSNE = pandas.read_csv(os.path.join(args.output, f"pyyamb.perplexity_{args.perplexity}.csv"))
		clusters = set(dfTSNE['cluster'])
		logger.info("Writing %s bins", len(clusters))
		frag_records = list(SeqIO.parse(args.assembly, "fasta"))
		for cluster in clusters:
			frag_names = list(dfTSNE[dfTSNE['cluster'] == cluster]['fragment'])
			output_bin = os.path.join(args.output, f"pyyamb.bin.{cluster}.fna")
			records = (x for x in frag_records if x.id in frag_names)
			write_records_to_fasta(records, output_bin, glue=True)

		logger.info("%s bins had been written", len(clusters))

	'''CheckM'''

	logger.info("Analysis finished")


if __name__ == '__main__':
	main()
