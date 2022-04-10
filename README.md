## DeepPerVar

By leveraging paired whole genome sequencing data and epigenetic functional assays in a population cohort, a  DeepPerVar is a multi-modal deep learning framework to predict genome-wide quantitative epigenetic signals and evaluate the functional consequence of noncoding variants on an individual level by quantifying their allelic difference on the prediction. By applying DeepPerVar to the ROSMAP cohort studying Alzheimer’s disease (AD), the web server can accurately predict genome-wide H3K9ac signals and DNA methylation ratio given DNA genomic sequence under reference and alternative alleles, and use the allelic difference as the score to evaluate the functional consequence of genetic variants associated with Alzheimer’s disease in a personal genome.

<center>

<div align=center><img width="800" height="500" src="https://raw.githubusercontent.com/alfredyewang/DeepPerVar/main/doc/DeepPerVar.jpeg"/></div>
</center>  


## DeepPerVar Webserver

We implement a webserver to predict genome-wide H3K9ac signals and DNA methylation ratio and the mutation effect on these two epigenetics signals. The webserver can be accessed from [link](http://35.202.146.70/). <br />
<br />

## Requirements and Installation

DeepPerVar is implemented by Python3.

- Python 3.8
- hdf5 == 1.10.4
- numpy >= 1.18.5
- pytorch ==1.7.1
- biopython=1.19.2

Install Samtools 1.15.1 follow the (instruction)[http://www.htslib.org/download/].

Download [Reference Genome (hg19)](https://drive.google.com/file/d/1X5PdUzaSVMKAkzysrv9kVtBbU3VGCkgF/view?usp=sharing), [DeepPerVar Models](https://drive.google.com/file/d/1Q_EzL_R4MLHSPYXKIqGUeXkDNx1yJ4WB/view?usp=sharing)

```
unzip Models.zip References.zip
```

Download DeepPerVar:
```
git clone https://github.com/alfredyewang/DeepPerVar
```
Install requirements, samtools can be downloaded and installed from [link](http://www.htslib.org/download/).
```
pip3 install -r requirements --user
```
## Usage
You can see the input arguments for DeepPerVar by help option:

```
usage: DeepPerVar.py [-h] [--prediction] [--epigenomics EPIGENOMICS] [--bed BED] [--model_dir <data_directory>] [--res_dir <data_directory>]

DeepPerVar: a multimodal deep learning framework for functional interpretation of genetic variants in personal genome

optional arguments:
  -h, --help            show this help message and exit
  --prediction          Use this option for predict DeepPerVar score
  --epigenomics EPIGENOMICS
                        Epigenetics, can be H3K9 or DNA_methylation
  --bed BED             The Bed file for predicts epigenetics and mutation effects
  --model_dir <data_directory>
                        The model directory for DeepPerVar
  --res_dir <data_directory>
                        The data directory for save results
```

### Input File Format
DeepPerVar takes UCSC Genome Browser BED file. Each line has 5 tab separated fields. The BED fields are:

- The first column: Chromosome name (hg19).
- The second column: Position of SNPs (hg19).
- The third column: The strand information.
- The forth column: reference allele.
- The fith column: alternative allele.
