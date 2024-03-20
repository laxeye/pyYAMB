#!/usr/bin/env python3
import argparse
import logging
import os
import pandas
import re
import seaborn
import hdbscan
import pycoverm
from statistics import mode
from Bio import SeqIO
from matplotlib import pyplot
from sklearn.manifold import TSNE
from pyyamb.cut_contigs import get_fragments
from pyyamb.map import map_reads
from pyyamb.map import view_mapping_file
from pyyamb.map import sort_mapping_file
from pyyamb.tetra import kmer_freq_table
from pyyamb.utils import write_records_to_fasta
from pyyamb.utils import check_files
from pyyamb import __version__


def parse_args():
	# logger = logging.getLogger()
	parser = argparse.ArgumentParser(prog="pyyamb",
		description=f"pyYAMB metagenome binner ver. {__version__}",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	general_args = parser.add_argument_group(title="General options", description="")
	general_args.add_argument("--task", required=True, default='all',
		choices=["all", "cut", "tetra", "map", "clustering", "make_bins"],
		help="Task of pipeline: cut (discard short contigs and cut longer), "
		+ "tetra (calculate frequencies of 4mers), "
		+ "map (map reads or process mapping file), "
		+ "clustering (perform tSNE and HDBSCAN to make clusters), "
		+ "make_bins (create FASTA files from fragments and clustering file).")
	general_args.add_argument("-o", "--output", type=str, required=True,
		help="Output directory.")
	general_args.add_argument("--min-length", default=1000, type=int,
		help="Minimum contig length.")
	general_args.add_argument("--fragment-length", default=10000, type=int,
		help="Target length of contig fragments.")
	general_args.add_argument("--perplexity", default=50, type=int,
		help="Perplexity parameter for tSNE.")
	general_args.add_argument("--min-cluster-size", "--mcs", default=30, type=int,
		help="HDBSCAN min_cluster_size parameter.")
	general_args.add_argument("--min-samples", "--ms", default=15, type=int,
		help="HDBSCAN min_samples parameter.")
	general_args.add_argument("--k-len", default=4, type=int,
		help="Length of k-mer to calculate their frequencies.")
	general_args.add_argument("-t", "--threads", type=int, default=1,
		help="Number of CPU threads to use (where possible).")
	general_args.add_argument('-u', '--write-unbinned', action='store_true',
		help="Write unbinned contigs to file.")
	general_args.add_argument("-d", "--debug", action='store_true',
		help="Print debug messages.")
	general_args.add_argument("--majority", action='store_true',
		help="All fragments from one contig go to the most represented bin.")
	general_args.add_argument("--force", action='store_true',
		help="Overwrite content of output directory.")
	general_args.add_argument('--version', action='version',
		version='%(prog)s {version}'.format(version=__version__))

	input_args = parser.add_argument_group(title="Input files and options")
	input_args.add_argument("-1", "--pe-1", nargs='*',
		help="First (left) paired-end reads, FASTQ [gzipped]. "
		+ "Space-separated if multiple.")
	input_args.add_argument("-2", "--pe-2", nargs='*',
		help="Second (right) paired-end reads, FASTQ [gzipped]. "
		+ "Space-separated if multiple.")
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
	input_args.add_argument("--clustered-data",
		help="Previously clustered data, only for \"make_bins\" task.")

	args = parser.parse_args()

	args.output = os.path.abspath(args.output)

	if not os.path.isdir(args.output):
		try:
			os.makedirs(args.output, exist_ok=args.force)
		except Exception as e:
			logging.error("Failed to create %s", args.output)
			raise e

	'''Check input'''
	for x in (
			args.assembly, args.kmers_data,
			args.coverage_data, args.clustered_data):
		if x:
			check_files([x])

	if args.single_end:
		args.single_end = list(map(os.path.abspath, args.single_end))
		check_files(args.single_end)

	if args.pe_1:
		args.pe_1 = list(map(os.path.abspath, args.pe_1))
		args.pe_2 = list(map(os.path.abspath, args.pe_2))
		check_files(args.pe_1 + args.pe_2)

	if args.mapping_file:
		args.mapping_file = list(map(os.path.abspath, args.mapping_file))
		check_files(args.mapping_file)

	return args


def create_logger(args):
	'''Create logger'''
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	if not args.debug:
		ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	fh = logging.FileHandler(os.path.join(args.output, "pyyamb.log"))
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	return logger


def make_nice_bam(args):
	'''Map reads and produce a sorted BAM file(s)'''
	sam_files = map_reads(args)
	bam_files = [view_mapping_file(args, sam, compress=False) for sam in sam_files]
	return [sort_mapping_file(args, bam_file) for bam_file in bam_files]


def main():
	args = parse_args()
	logger = create_logger(args)
	mg_data = pandas.DataFrame()
	logger.info("Analysis started")

	'''Cut contigs'''
	if args.task in ("cut", "all"):
		try:
			fragments = get_fragments(
				args.assembly,
				args.fragment_length,
				args.min_length)
			args.assembly = write_records_to_fasta(
				fragments,
				os.path.join(args.output, 'fragments.fasta'))
			logger.info("Contigs fragmented.")
		except Exception as e:
			logger.error("Error during contig fragmentation.")
			raise e

	'''Count kmers'''
	if args.task in ("tetra", "all"):
		try:
			logger.info("Counting kmers")
			k_data = kmer_freq_table(args.assembly, args.k_len, args.threads)
			kmers_data_file = os.path.join(args.output, "kmers.csv")
			k_data.to_csv(kmers_data_file)
			logger.info('Wrote k-mer frequencies to %s.', kmers_data_file)
		except Exception as e:
			logger.error("Error during kmer frequency calculation.")
			raise e

	'''Map reads with minimap2'''
	if (args.task == "map" or
		(args.task == 'all' and (args.pe_1 or args.single_end or args.mapping_file))):
		if args.mapping_file:
			sorted_bam_files = args.mapping_file
		else:
			sorted_bam_files = make_nice_bam(args)
		logger.info("Extracting coverage")
		samples = [os.path.splitext(os.path.basename(x))[0] for x in sorted_bam_files]
		contigs, mean_covs = pycoverm.get_coverages_from_bam(
			sorted_bam_files, contig_end_exclusion=75, min_identity=0.97,
			trim_lower=0.0, trim_upper=0.0, threads=args.threads)
		cov_data = pandas.DataFrame(mean_covs, index=contigs, columns=samples)
		cov_data.index.rename('fragment', inplace=True)
		cov_data.to_csv(os.path.join(args.output, "coverage.csv"))
		logger.info("Processing of mapping file finished")

	'''Cluster data with HDBSCAN'''
	if args.task == "clustering":
		'''Read k-mers from disk. Frequencies, not z-scores'''
		if args.kmers_data and os.path.isfile(args.kmers_data):
			k_data = pandas.read_csv(args.kmers_data, index_col=0)
			logger.info('Read k-mer frequencies from %s.', args.kmers_data)
		else:
			logger.error("File with k-mers not found")
			raise FileNotFoundError(args.kmers_data)

		'''Read coverage depth of fragments from disk. Not z-scores'''
		if args.coverage_data and os.path.isfile(args.coverage_data):
			cov_data = pandas.read_csv(args.coverage_data)
			logger.info('Read fragment coverage depth from %s.', args.coverage_data)
		else:
			logger.error("File fragment coverage depth not found")
			raise FileNotFoundError(args.coverage_data)

	if args.task in ("clustering", "all"):
		'''Merge data, produce Z-scores and dump to file'''
		mg_data = k_data.merge(cov_data, how='left', on='fragment')
		z_scores = mg_data.iloc[:, :2].join(
			mg_data.iloc[:, 2:].apply(
				lambda x: (x - x.mean()) / x.std(),
				axis=0
			)
		)
		z_scores.to_csv(os.path.join(args.output, "z-scored_data.csv"))

		logger.info('tSNE data reduction.')
		tsne_pca = TSNE(
			init='pca',
			perplexity=args.perplexity).fit_transform(z_scores.iloc[:, 2:])
		dfTSNE = pandas.DataFrame.join(
			pandas.DataFrame(tsne_pca, columns=['tsne1', 'tsne2'], index=k_data.index),
			z_scores[['fragment', 'length']]
		)

		logger.info('HDBSCAN data clustering.')
		clusterer = hdbscan.HDBSCAN(
			min_cluster_size=args.min_cluster_size,
			min_samples=args.min_samples
		)
		cluster_labels = clusterer.fit_predict(dfTSNE[['tsne1', 'tsne2']])
		dfTSNE = dfTSNE.join(
			pandas.DataFrame(cluster_labels, columns=['cluster'])
		)
		fname_prefix = "_".join([
			"pyyamb.tp", str(args.perplexity),
			"mcl", str(args.min_cluster_size),
			"ms", str(args.min_samples)])
		dfTSNE.to_csv(os.path.join(args.output, f"{fname_prefix}.csv"))
		clusters = set(dfTSNE['cluster'])
		logger.info("HDBSCAN found %s clusters", len(clusters))

		logger.info("Producing graphics")
		pal = seaborn.color_palette(palette="Set3", n_colors=len(clusters))
		seaborn.relplot(
			data=dfTSNE, x='tsne1', y='tsne2', size='length',
			alpha=0.1, hue=dfTSNE["cluster"].astype("category"),
			palette=pal
		)
		pyplot.savefig(os.path.join(args.output, f"{fname_prefix}.png"), dpi=300)
		pyplot.savefig(os.path.join(args.output, f"{fname_prefix}.svg"))

	if args.task in ["make_bins", "clustering", "all"]:
		if args.task == "make_bins":
			dfTSNE = pandas.read_csv(os.path.join(args.clustered_data))
		frag_records = list(SeqIO.parse(args.assembly, "fasta"))
		bin_dir = os.path.join(args.output, "bins")
		os.makedirs(bin_dir, exist_ok=args.force)

		if args.majority:
			pattern = re.compile(r'^(.+)_frag_(\d+)$')
			buffer = []
			c_prev = None
			bin_N = []
			out_dict = dict()
			with open(os.path.join(args.output, "pyyamb.contig2bin.txt"), 'w') as c2b:
				for i in dfTSNE.index:
					row = dfTSNE.iloc[i]
					cluster, frag = row['cluster'], row['fragment']
					m = pattern.match(frag)
					if m:
						contig = m.group(1)
						if contig == c_prev:
							buffer.append(frag)
							bin_N.append(cluster)
						else:
							if len(buffer) > 0:
								final_bin = mode(bin_N)
								out_dict[final_bin] = out_dict.get(final_bin, []) + buffer
								c2b.write(f'{c_prev},{final_bin}\n')
							buffer = [frag]
							bin_N = [cluster]
							c_prev = contig
					else:
						if len(buffer) > 0:
							final_bin = mode(bin_N)
							out_dict[final_bin] = out_dict.get(final_bin, []) + buffer
							c2b.write(f'{c_prev},{final_bin}\n')
						else:
							out_dict[cluster] = out_dict.get(cluster, []) + [frag]
							c2b.write(f'{frag},{final_bin}\n')
						buffer = []
			bin_count = 0
			for k, frag_names in out_dict.items():
				if int(k) != -1:
					bin_file = os.path.join(bin_dir, f"pyyamb.bin.{k}.fna")
				elif args.write_unbinned and int(k) == -1:
					bin_file = os.path.join(bin_dir, "pyyamb.unbinned.fna")
				else:
					continue
				records = (x for x in frag_records if x.id in frag_names)
				write_records_to_fasta(records, bin_file, glue=True)
				bin_count += 1
			logger.info("%s bins had been written", bin_count)

		else:
			clusters = set(dfTSNE['cluster'])
			if not args.write_unbinned and -1 in clusters:
				clusters.remove(-1)
			logger.info("Writing %s bins", len(clusters))
			for cluster in clusters:
				frag_names = list(dfTSNE[dfTSNE['cluster'] == cluster]['fragment'])
				if args.write_unbinned and cluster == -1:
					output_bin = os.path.join(bin_dir, "pyyamb.unbinned.fna")
				else:
					output_bin = os.path.join(bin_dir, f"pyyamb.bin.{cluster}.fna")
				records = (x for x in frag_records if x.id in frag_names)
				write_records_to_fasta(records, output_bin, glue=True)

			logger.info("%s bins had been written", len(clusters))

	logger.info("Analysis finished")


if __name__ == '__main__':
	main()
