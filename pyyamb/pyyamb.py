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
		choices=["all", "assemble", "cut", "tetra", "tsne", "map", "bin", "checkm"])
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
	input_args.add_argument("-s", "--single-end", nargs='+',
		help="Sinle-end reads, FASTQ. Space-separated if multiple.")
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
	parser.add_argument("--perplexity", default=50, type=int,
		help="Perplexity parameter for tSNE.")
	parser.add_argument("--kmers-data",
		help="Previously calculated kmer-freqs.")

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
		if os.path.isfile(args.assembly):
			args.assembly = os.path.abspath(args.assembly)
		else:
			raise FileNotFoundError(args.assembly)

	if args.single_end:
		args.single_end = list(map(os.path.abspath, args.single_end))
	if args.pe_1:
		args.pe_1 = list(map(os.path.abspath, args.pe_1))
		args.pe_2 = list(map(os.path.abspath, args.pe_2))
		for f in args.pe_1 + args.pe_2:
			if not os.path.isfile(f):
				raise FileNotFoundError(f"Input file {f} not found!")

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
	logger = createLogger()
	args = parse_args()
	logger.info("Analysis started")

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

	if args.task == "map" or args.task == 'all':
		'''Map reads with minimap2'''
		target = os.path.join(args.output, "mapping.sam")
		mapping_bam_file = os.path.join(args.output, 'mapping.bam')
		if args.single_end:
			# Only first single-end reads file is currently mapped
			reads = [args.single_end[0]]
		if args.pe_1 and args.pe_1:
			# Only one pair of reads is currently mapped
			reads = list(zip(args.pe_1, args.pe_2))[0]
		mapping_sam_file = map_reads(args.assembly, reads, target, args.threads)
		logger.info("Converting mapping file")
		with open(mapping_bam_file, 'wb') as h:
			h.write(pysam.view('-@', str(args.threads), '-u', '-F', '4', mapping_sam_file))
		sorted_mapping_bam_file = os.path.join(args.output, 'mapping.sorted.bam')
		logger.info("Sorting mapping file")
		pysam.sort('-@', str(args.threads), '-o', sorted_mapping_bam_file, mapping_bam_file)
		logger.info("Indexing mapping file")
		pysam.samtools.index(sorted_mapping_bam_file)
		logger.info("Extracting coverage")
		bam_handle = pysam.AlignmentFile(sorted_mapping_bam_file)
		contig_list = [record.id for record in SeqIO.parse(args.assembly, 'fasta')]
		d = {contig: mean(bam_handle.count_coverage(contig=contig)) for contig in contig_list}
		bam_handle.close()
		coverage_df = pandas.DataFrame.from_dict(d, orient='index', columns=['coverage'])
		coverage_df.to_csv(os.path.join(args.output, "coverage.csv"))
		logger.info("Processing of mapping file finished")

	'''Count kmers'''
	if args.task == "tetra" or args.task == "all":
		try:
			if args.kmers_data and os.path.isfile(args.kmers_data):
				freqs_df = pandas.read_csv(args.kmers_data, index_col=0)
				logger.info('Read kmer frequencies from %s.', args.kmers_data)
			else:
				logger.info("Counting kmers")
				freqs_df = kmer_freq_table(args.assembly)
				args.kmers_data = os.path.join(args.output, "kmers.csv")
				freqs_df.to_csv(args.kmers_data)
				logger.info('Kmer frequencies written to %s.', args.kmers_data)

		except Exception as e:
			logger.error("Error during kmer frequency calculation.")
			raise e

	if args.task == "tsne" or args.task == "all":
		if args.kmers_data and os.path.isfile(args.kmers_data):
			freqs_df = pandas.read_csv(args.kmers_data, index_col=0)
			logger.info('Read kmer frequencies from %s.', args.kmers_data)
		else:
			logger.error("File with k-mers not found")
			raise FileNotFoundError(args.kmers_data)

		try:
			freqs_df.join(coverage_df, how='inner')
		except Exception as e:
			logger.debug("Coverage missing, skipping.")

		logger.info('tSNE analysis started.')
		tsne_pca = TSNE(init='pca', perplexity=args.perplexity).fit_transform(freqs_df.iloc[:, 1:])
		logger.info('tSNE analysis complete.')
		dfTSNE = pandas.DataFrame.join(
			pandas.DataFrame(tsne_pca, columns=['tsne1', 'tsne2'], index=freqs_df.index),
			freqs_df['length']
		)

		logger.info('HDBSCAN analysis started.')
		clusterer = hdbscan.HDBSCAN(min_cluster_size=25)
		cluster_labels = clusterer.fit_predict(dfTSNE[['tsne1', 'tsne2']])
		logger.info('HDBSCAN analysis completed.')
		dfTSNE = dfTSNE.join(pandas.DataFrame(
			cluster_labels, index=freqs_df.index, columns=['cluster']
			))
		dfTSNE.to_csv(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.csv"))
		clusters = set(dfTSNE['cluster'])
		logger.info("%s clusters found", len(clusters))

		frag_records = list(SeqIO.parse(args.assembly, "fasta"))

		for cluster in clusters:
			dfCluster = dfTSNE[dfTSNE['cluster'] == cluster]
			dfCluster.to_csv(os.path.join(args.output, f"bin.{cluster}.csv"))
			output_bin = os.path.join(args.output, f"bin.{cluster}.fna")
			sequences = [x for x in frag_records if x.id in list(dfCluster.index)]
			write_records_to_fasta(sequences, output_bin)

		logger.info("%s bins had been written", len(clusters))

		pal = seaborn.color_palette(palette="Set2", n_colors=len(clusters))
		seaborn.relplot(
			data=dfTSNE, x='tsne1', y='tsne2', size='length',
			alpha=0.2, hue=dfTSNE["cluster"].astype("category"),
			palette=pal
		)
		pyplot.savefig(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.png"), dpi=300)
		pyplot.savefig(os.path.join(args.output, f"tetra.tsne.{args.perplexity}.svg"))

	'''Process data'''
	logger.info("Analysis finished")


if __name__ == '__main__':
	main()
