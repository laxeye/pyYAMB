# pyYAMB

[![Python versions](https://img.shields.io/pypi/pyversions/pyyamb.svg)](https://pypi.org/project/pyYAMB/)
[![PyPI version](https://img.shields.io/pypi/v/pyyamb.svg)](https://pypi.org/project/pyYAMB/)
[![Anaconda version](https://anaconda.org/laxeye/pyyamb/badges/version.svg)](https://anaconda.org/laxeye/pyyamb/)

pyYAMB is an implementation of [YAMB](https://github.com/laxeye/YAMB/) (Yet another metagenome binner) on Python (>=3.8). Originally YAMB was described in the preprint https://www.biorxiv.org/content/10.1101/521286.abstract and it's main point is the use of tSNE and HDBSCAN to process tetramer frequencies and coverage depth of metagenome fragments. pyYAMB strives for parallel computing wherever possible.

### pyYAMB data processing includes

* contig filtering and fragmentation
* read mapping with [minimap2](https://github.com/lh3/minimap2)
* mapping files processing with pysam and coverage depth extraction with [pycoverm](https://github.com/apcamargo/pycoverm)
* k-mer (by default tetramer) frequency calculation
* data dimensions reduction with tSNE
* data clustering with HDBSCAN
* writing bins to FASTA
* writing plots to PNG and SVG

#### Possible features in *far* future

* read processing
* metagenome assembly
* bin QC

## How to start

### Installation

#### PyPI

pyYAMB is available at PyPI and may be installed with:

`pip install pyYAMB`

Also yo need to install dependencies (see below).

#### GitHub

Another way (*not recommended*) is to clone the repository

`git clone https://github.com/laxeye/pyYAMB.git` or `gh repo clone laxeye/pyYAMB`

and run

`python setup.py install` or `pip install .`

It installs pyYAMB and python libraries. Problems may appear with *hdbscan* module and *cython*. Just reinstall *hdbscan* using `pip install hdbscan` and try again `python setup.py install`.

#### Dependencies

If you installed pyYAMB from PyPI or GitHub, you need to install dependencies: minimap2 and samtools (e.g. using conda).

`conda install -c bioconda minimap2 "samtools>=1.9"`

#### Conda

Currently only outdated versions of pyYAMB are available at Anaconda. They may be installed with all dependencies:

`conda install -c laxeye pyyamb` or `mamba install -c laxeye pyyamb`

### Usage

pyYAMB entry point is the all-in-one command `pyyamb`. pyYAMB has two dozens of arguments, their description is available after running `pyyamb -h`

You may start from metagenome assembly and processed (quality trimmed etc.) reads, e.g.:

`pyyamb --task all -1 Sample_1.R1.fastq.gz Sample_2.R1.fastq.gz -2 Sample_1.R2.fastq.gz Sample_2.R2.fastq.gz -i assembly.fasta -o results/will/be/here --threads 8`

After completion bins could be found in *bins* subfolder in output folder. "-1" bin collects unbinned sequences. Quality check of resulting bins is strongly recommended, You may use [CheckM](https://github.com/Ecogenomics/CheckM) or [CheckM2](https://github.com/chklovski/CheckM2).

## Results and benchmarks

pyYAMB was tested on low complexity data set for the 1st CAMI challenge (simulated Illumina HiSeq data, small insert size).

- Number of samples: 1
- Total Size: 15 Gbp
- Read length: 2x150 bp
- Insert size mean: 270 bp
- Insert size stddev: 27 bp

The run took 12 minutes and 17 seconds on AMD Ryzen 3900X using 8 threads. Completeness and purity results are given below:

|Property|Value|
|--------|-----|
|Average completeness (bp)|94.4%|
|Average completeness (seq)|84.2%|
|Average purity (bp)|67.1%|
|Average purity (seq)|56.5%|


Earlier YAMB showed quality comparable with [CONCOCT](https://github.com/BinPro/CONCOCT) binner (see the [preprint](https://www.biorxiv.org/content/10.1101/521286.abstract) for details).


## References

Van Der Maaten, L. (2014). Accelerating t-SNE using tree-based algorithms. The Journal of Machine Learning Research, 15(1), 3221-3245.

Campello, R. J., Moulavi, D., & Sander, J. (2013, April). Density-based clustering based on hierarchical density estimates. In Pacific-Asia conference on knowledge discovery and data mining (pp. 160-172). Springer, Berlin, Heidelberg.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. https://dx.doi.org/10.1093/bioinformatics/bty191

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055. https://dx.doi.org/10.1101/gr.186072.114

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. https://doi.org/10.1093/bioinformatics/btp163
