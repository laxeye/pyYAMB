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
from sklearn.manifold import TSNE
from pyyamb.cut_contigs import get_fragments
from pyyamb.map import map_reads
from pyyamb.utils import write_records_to_fasta
from pyyamb.tetra import kmer_freq_table
from numpy import mean


def parse_args():
	logger = logging.getLogger()
	parser = argparse.ArgumentParser(prog="pyYAMB",
		description="pyYAMB metagenome binner")

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", help="Task of pipeline", required=True,
		choices=["all", "cut", "tetra", "clustering", "map"])
	general_args.add_argument("--mapper", help="Mapping software, currently minimap2 only",
		choices=["minimap2"], default="minimap2")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible), default 1.")
	general_args.add_argument("-m", "--memory-limit", type=int, default=16,
		help="Memory limit in GB, default 16.")
	general_args.add_argument("-o", "--output", type=str, required=True,
		help="Output directory")
	general_args.add_argument("--min-length", default=1000, type=int,
		help="minimum contig length")
	general_args.add_argument("--fragment-length", default=10000, type=int,
		help="Target length of contig fragments")
	general_args.add_argument("--perplexity", default=50, type=int,
		help="Perplexity parameter for tSNE.")
	general_args.add_argument("--k-len", default=4, type=int,
		help="Length of k-mer, default: 4.")

	input_args = parser.add_argument_group(title="Input files and options")
	input_args.add_argument("-1", "--pe-1", nargs='+',
		help="First (left) paired-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='+',
		help="Second (right) paired-end reads, FASTQ. Space-separated if multiple.")
	input_args.add_argument("-s", "--single-end", nargs='+',
		help="Sinle-end reads, FASTQ. Space-separated if multiple.")
	input_args.add_argument("-i", "--assembly",
		help="Previously assembled metagenome.")
	input_args.add_argument("--kmers-data",
		help="Previously calculated kmer-freqs.")

	'''asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler",
		default="unicycler", choices=["spades", "megahit"],
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
	if args.assembly:
		args.assembly = os.path.abspath(args.assembly)
		if not os.path.isfile(args.assembly):
			raise FileNotFoundError(args.assembly)

	if args.single_end:
		args.single_end = list(map(os.path.abspath, args.single_end))

	if args.pe_1:
		args.pe_1 = list(map(os.path.abspath, args.pe_1))
		args.pe_2 = list(map(os.path.abspath, args.pe_2))
		for f in args.pe_1 + args.pe_2:
			if not os.path.isfile(f):
				raise FileNotFoundError(f"Input file {f} not found!")

	if args.kmers_data:
		args.kmers_data = os.path.abspath(args.kmers_data)
		if not os.path.isfile(args.kmers_data):
			raise FileNotFoundError(f"Input file {f} not found!")

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
	with pysam.AlignmentFile(bamfile) as bam_handle:
		contig_list = [record.id for record in SeqIO.parse(fasta, 'fasta')]
		d = {contig: mean(bam_handle.count_coverage(contig=contig)) for contig in contig_list}
	return pandas.DataFrame.from_dict(d, orient='index', columns=['coverage'])


def main():
	logger = create_logger()
	args = parse_args()
	logger.info("Analysis started")
	mg_data = pandas.DataFrame()

	'''Process reads'''

	'''Assemble with spades or megahit'''

	if args.task == "cut" or args.task == "all":
		try:
			fragments = get_fragments(args.assembly, args.fragment_length, args.min_length)
			args.assembly = write_records_to_fasta(fragments, os.path.join(args.output, 'fragments.fna'))
			logger.info("Contigs fragmented.")
		except Exception as e:
			logger.error("Error during contig fragmentation.")
			raise e

	'''Count kmers'''
	if args.task == "tetra" or args.task == "all":
		if args.kmers_data and os.path.isfile(args.kmers_data):
			mg_data = pandas.read_csv(args.kmers_data, index_col=0)
			logger.info('Read k-mer frequencies from %s.', args.kmers_data)
		else:
			try:
				logger.info("Counting kmers")
				mg_data = kmer_freq_table(args.assembly, args.k_len)
				args.kmers_data = os.path.join(args.output, "kmers.csv")
				mg_data.to_csv(args.kmers_data)
				logger.info('Wrote k-mer frequencies to %s.', args.kmers_data)
			except Exception as e:
				logger.error("Error during kmer frequency calculation.")
				raise e

	'''Map reads with minimap2'''
	if args.task == "map" or (args.task == 'all' and (args.pe_1 or args.single_end)):
		mapping_sam_file = map_reads(args)

		logger.info("Converting mapping file")
		mapping_bam_file = os.path.join(args.output, 'mapping.bam')
		with open(mapping_bam_file, 'wb') as h:
			h.write(pysam.view('-@', str(args.threads), '-u', '-F', '4', mapping_sam_file))
		os.remove(mapping_sam_file)

		logger.info("Sorting mapping file")
		sorted_mapping_bam_file = os.path.join(args.output, 'mapping.sorted.bam')
		pysam.sort('-@', str(args.threads), '-o', sorted_mapping_bam_file, mapping_bam_file)
		os.remove(mapping_bam_file)

		logger.info("Indexing mapping file")
		pysam.samtools.index(sorted_mapping_bam_file)

		logger.info("Extracting coverage")
		cov_data = extract_coverage(sorted_mapping_bam_file, args.assembly)
		cov_data.to_csv(os.path.join(args.output, "coverage.csv"))
		mg_data = mg_data.join(cov_data)
		logger.info("Processing of mapping file finished")

	'''Cluster data with HDBSCAN'''
	if args.task == "clustering" or args.task == "all":
		'''!Remove excess read from disk!'''
		if args.kmers_data and os.path.isfile(args.kmers_data):
			freqs_df = pandas.read_csv(args.kmers_data, index_col=0)
			logger.info('Read k-mer frequencies from %s.', args.kmers_data)
		else:
			logger.error("File with k-mers not found")
			raise FileNotFoundError(args.kmers_data)

		logger.info('tSNE data reduction.')
		tsne_pca = TSNE(init='pca', perplexity=args.perplexity).fit_transform(mg_data.iloc[:, 1:])
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
		dfTSNE.to_csv(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.csv"))
		clusters = set(dfTSNE['cluster'])
		logger.info("HDBSCAN found %s clusters", len(clusters))

		frag_records = list(SeqIO.parse(args.assembly, "fasta"))

		for cluster in clusters:
			dfCluster = dfTSNE[dfTSNE['cluster'] == cluster]
			'''dfCluster.to_csv(os.path.join(args.output, f"bin.{cluster}.csv"))'''
			output_bin = os.path.join(args.output, f"bin.{cluster}.fna")
			sequences = [x for x in frag_records if x.id in list(dfCluster.index)]
			write_records_to_fasta(sequences, output_bin)

		logger.info("%s bins had been written", len(clusters))

		pal = seaborn.color_palette(palette="Set3", n_colors=len(clusters))
		seaborn.relplot(
			data=dfTSNE, x='tsne1', y='tsne2', size='length',
			alpha=0.2, hue=dfTSNE["cluster"].astype("category"),
			palette=pal
		)
		pyplot.savefig(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.png"), dpi=300)
		pyplot.savefig(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.svg"))

	'''CheckM'''
	
	logger.info("Analysis finished")


if __name__ == '__main__':
	main()
