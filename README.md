# MQ2Gpipe
 A simple Snakemake workflow for genomic variants calling and annotation form Raw NGS data of multiple individuals.

Author: Peng Liu

Email: sxliulian2012@hotmail.com

# Installation

**Dependencies**

* conda > v3.8.5
* snakemake > v5.7.0
* fastp (version: 0.20.1)
* bwa (version: 2.2.1)
* samtools (version: 1.9)
* gatk > v4.0
* ...

Conda can be downloaded as part of the Anaconda or the Miniconda plattforms. We recommend to install miniconda3. Using Linux you can get it with:

```sh
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```


Snakemake can be install with:

```sh
conda create -c conda-forge -c bioconda -n MQ2Gpipe_env python=3 snakemake
```

Detail installation guide from: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

other packages can be installed by conda:

```sh
# install normal softwares
conda install -n MQ2Gpipe_env -c bioconda -c conda-forge fastp bwa samtools bedtools vcfotools bcftools plink beagle gatk4 ensembl-vep 
```

MQ2Gpipe.smk 

# Usage

## Data preparation

```sh
# 1. download maize reference genome fasta
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
wget -c ftp://ftp.gramene.org/pub/gramene/release-63/gtf/zea_mays/Zea_mays.B73_RefGen_v4.48.gtf.gz
# 2. uncompression
gzip -cd Zea_mays.AGPv4.dna.toplevel.fa.gz |cut -f1 > /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4/zma_dna_v4.fa
gzip -cd Zea_mays.B73_RefGen_v4.48.gtf.gz > /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4/zma.v4.48.gtf

## build bwa index
bwa index genome.fa
## ctrate genome dict
samtools dict -o zma_dna_v4.fa.dict zma_dna_v4.fa

```
## Configure

**Example:**
```yaml
### Path to an uncompressed FASTA file with all transcript sequences.
# fasta
dna: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/zma_v4.dna.fa
# annotation
gtf: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/Zea_mays.B73_RefGen_v4.50.gtf
# reference vcf
vcf: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/Variation/zea_mays_v4.vcf.gz
tbi: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/Variation/zea_mays_v4.vcf.gz.tbi
vep_db: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/VEP_Cache/zea_mays.50.vep_db
vep_ver: 50
# index of genome
genome_dict: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/zma_v4.dna.fa.dict
genome_bwa_index: /MaizeLab/auhpc1/DataBase/Species/Zea_Mays/B73v4.50/BWA_Index/zma_v4.dna.fa

## Path to a folder where intermediate files will be written.
output_dir: /MaizeLab/auhpc1/Project/MQ2G_Output/

# Path to a YAML file with samples and their corresponding FASTQ files.
sample_list: /MaizeLab/auhpc1/Project/sample_list.yaml

```

## Samples

**Example:**

```yaml
# sample.yaml
SampleA:
- /Path/to/sampleA_R1.fastq.gz
- /Path/to/sampleA_R2.fastq.gz
SampleB:
- /Path/to/sampleB_R1.fastq.gz
- /Path/to/sampleB_R2.fastq.gz
SampleC:
- /Path/to/sampleC_R1.fastq.gz
- /Path/to/sampleC_R2.fastq.gz

```


## Run pipeline

```sh
source activate MQ2Gpie_env
snakemake -p -s MQ2Gpie.smk -j <threads> --latency-wait 20 
```

* threads: The threads number for use. Less than max threads.

# FAQs

**Q:** error while loading shared libraries: libreadline.so.6

**A:** link `libreadline.so.?` to `libreadline.so.6`. example: `ln -s  /opt/miniconda3/envs/feelnc_env/lib/libreadline.so.8  /opt/miniconda3/envs/feelnc_env/lib/libreadline.so.6`

# Further

More features will be added soon ...