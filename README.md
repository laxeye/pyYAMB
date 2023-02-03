# pyYAMB

[![Python versions](https://img.shields.io/pypi/pyversions/pyyamb.svg)](https://pypi.org/project/pyYAMB/)
[![PyPI version](https://img.shields.io/pypi/v/pyyamb.svg)](https://pypi.org/project/pyYAMB/)

pyYAMB is an implementation of [YAMB](https://github.com/laxeye/YAMB/) (Yet another metagenome binner) on Python (>=3.8). Originally YAMB was described in the preprint https://www.biorxiv.org/content/10.1101/521286.abstract and it's main point is the use of tSNE and HDBSCAN to process tetramer frequencies and coverage depth of metagenome fragments. pyYAMB strives for parallel computing wherever possible, currently coverage depth extraction is single threaded and takes the most time.

### pyYAMB data processing includes

* contig filtering and fragmentation
* read mapping with minimap2
* mapping files processing and coverage depth extraction with pysam
* k-mer (by default tetramer) frequency calculation
* data dimensions reduction with tSNE
* data clustering with HDBSCAN
* writing bins to FASTA

#### Possible features in *far* future

* read processing
* metagenome assembly
* bin QC

## How to start

**Warning!** *pyYAMB now is in alpha-testing and may be unstable, use it at Your own risk.*

### Installation

#### PyPI

pyYAMB is available at PyPI and may be installed with:

`pip install pyYAMB`

#### GitHub

Another way is to clone the repository

`git clone https://github.com/laxeye/pyYAMB.git` or `gh repo clone laxeye/pyYAMB`

and run

`python setup.py install`

It installs pyYAMB and python libraries. Problems may appear with *hdbscan* module and *cython*. Just reinstall *hdbscan* using `pip install hdbscan` and try again `python setup.py install`.

#### Dependencies

Also you need to install dependencies: minimap2. Samtools should be installed automatically during pysam installation, otherwise please install it (e.g. using conda). 

Conda package will be available soon.

### Usage

pyYAMB entry point is the all-in-one command `pyyamb`. pyYAMB has two dozens of arguments, their description is available after running `pyyamb -h`

You may start from metagenome assembly and processed (quality trimmed etc.) reads, e.g.:

`pyyamb --task all -1 Sample_1.R1.fastq.gz Sample_2.R1.fastq.gz -2 Sample_1.R2.fastq.gz Sample_2.R2.fastq.gz -i assembly.fasta -o results/will/be/here --threads 8`

After completion bins could be found in *bins* subfolder in output folder. "-1" bin collects unbinned sequences.

## Results and benchmarks

pyYAMB will be tested on CAMI dataset soon. YAMB showed quality compared with [CONCOCT](https://github.com/BinPro/CONCOCT) binner (see the [preprint](https://www.biorxiv.org/content/10.1101/521286.abstract) for details).


## References

Van Der Maaten, L. (2014). Accelerating t-SNE using tree-based algorithms. The Journal of Machine Learning Research, 15(1), 3221-3245.

Campello, R. J., Moulavi, D., & Sander, J. (2013, April). Density-based clustering based on hierarchical density estimates. In Pacific-Asia conference on knowledge discovery and data mining (pp. 160-172). Springer, Berlin, Heidelberg.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. https://dx.doi.org/10.1093/bioinformatics/bty191

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055. https://dx.doi.org/10.1101/gr.186072.114

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. https://doi.org/10.1093/bioinformatics/btp163
