#!/usr/bin/env python3
import argparse
import logging
import os
import pandas
import seaborn
import hdbscan
from matplotlib import pyplot
from sklearn.manifold import TSNE
from pyyamb.cut_contigs import get_fragments
from pyyamb.map import map_reads
from pyyamb.utils import write_records_to_fasta
from pyyamb.tetra import kmer_freq_table


def parse_args():
	parser = argparse.ArgumentParser(prog="pyYAMB",
		description="pyYAMB metagenome binner")

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", help="Task of pipeline",
		choices=["assemble", "cut", "tetra", "map", "bin", "checkm"], required=True)
	general_args.add_argument("--mapper", help="Mapping software",
		choices=["minimap2", "bwa", "bowtie2"], default="minimap2")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible), default 1.")
	general_args.add_argument("-m", "--memory-limit", type=int, default=16,
		help="Memory limit in GB, default 16.")
	general_args.add_argument("-o", "--output", type=str, required=True,
		help="Output directory")

	input_args = parser.add_argument_group(title="Input files and options")
	input_args.add_argument("-1", "--pe-1", nargs='+',
		help="First (left) paired-end reads, FASTQ [gzipped]. Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='+',
		help="Second (right) paired-end reads, FASTQ. Space-separated if multiple.")
	input_args.add_argument("-i", "--assembly",
		help="Previously assembled metagenome.")

	asly_args = parser.add_argument_group(title="Assembly settings")
	asly_args.add_argument("-a", "--assembler",
		default="unicycler", choices=["spades", "megahit"],
		help="Assembler: 'spades' or 'megahit'.")
	asly_args.add_argument("--spades-correction", action="store_true",
		help="Perform short read correction by SPAdes (not recommended).")
	asly_args.add_argument("--spades-k-list",
		help="SPAdes: List of kmers, comma-separated even numbers e.g. '21,33,55,77'")

	parser.add_argument("--min-length", default=1000, type=int,
		help="minimum contig length")
	parser.add_argument("--fragment-length", default=10000, type=int,
		help="Target length of contig fragments")

	args = parser.parse_args()

	args.output = os.path.abspath(args.output)

	'''Check input'''

	if args.assembly:
		if os.path.isfile(args.assembly):
			args.assembly = os.path.abspath(args.assembly)
		else:
			raise FileNotFoundError(args.assembly)

	if args.pe_1:
		if os.path.isfile(args.pe_1) and os.path.isfile(args.pe_2):
			args.pe_1 = os.path.abspath(args.pe_1)
			args.pe_2 = os.path.abspath(args.pe_2)
		else:
			raise FileNotFoundError(f"Input files {args.pe_1} and {args.pe_2} not found!")

	return args


def createLogger():
	'''Init logger'''
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	return logger


def main():
	args = parse_args()

	logger = createLogger()

	'''Process reads'''

	'''Assemble with spades, megahit, ?'''

	if args.task == "cut":
		try:
			fragments = get_fragments(args.assembly, args.fragment_length, args.min_length)
			fragments_file = write_records_to_fasta(fragments, os.path.join(args.output, 'fragments.fna'))
			logger.info("Contigs fragmented.")
		except Exception as e:
			logger.error("Error during contig fragmentation.")
			raise e

	'''Count kmers'''
	if args.task == "tetra":
		try:
			perplexity = 50
			freqs_df = kmer_freq_table(args.assembly)
			freqs_file = os.path.join(args.output, "kmers.csv")
			freqs_df.to_csv(freqs_file)
			logger.info('Kmer frequencies written to %s.', freqs_file)
			tsne_pca = TSNE(init='pca', perplexity=perplexity).fit_transform(freqs_df.iloc[:, 1:])
			logger.info('tSNE analysis complete.')
			dfTSNE = pandas.DataFrame.join(
				pandas.DataFrame(tsne_pca, columns=['tsne1', 'tsne2'], index=freqs_df.index),
				freqs_df['length']
			)
			clusterer = hdbscan.HDBSCAN(min_cluster_size=25)
			cluster_labels = clusterer.fit_predict(dfTSNE[['tsne1', 'tsne2']])
			logger.info('HDBSCAN analysis complete.')
			dfTSNE = dfTSNE.join(pandas.DataFrame(
				cluster_labels, index=freqs_df.index, columns=['cluster']
				))
			dfTSNE.to_csv("tetra.tsne.csv")
			clusters = set(dfTSNE['cluster'])
			logger.info("%s clusters found", len(clusters))
			for cluster in clusters:
				dfTSNE[dfTSNE['cluster'] == cluster].to_csv(f"bin.{cluster}.csv")
			seaborn.relplot(data=dfTSNE, x='tsne1', y='tsne2', size='length',
				alpha=0.2, hue=dfTSNE["cluster"].astype("category"),
				palette=seaborn.color_palette(palette="Set2", n_colors=len(clusters))
			)
			pyplot.savefig(f"tetra.tsne.{perplexity}.png", dpi=300)
			pyplot.savefig(f"tetra.tsne.{perplexity}.svg")

		except Exception as e:
			logger.error("Error during kmer frequency calculation.")
			raise e

	'''Map with minimap2'''
	if args.task == "map":
		target = os.path.join(args.output, "mapping.sam")
		if map_reads(args.assembly, (args.pe_1, args.pe_2), target, args.threads):
			logger.info("Reads mapped.")
		else:
			logger.error("Error during mapping.")

	'''Extract coverage'''

	'''Process data'''


if __name__ == '__main__':
	main()
