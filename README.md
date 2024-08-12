# seq2squiggle-benchmark
Collection of scripts for benchmarking of seq2squiggle with other simulators and experimental data.

## Setup and running
This repo contains an automated workflow to benchmark nanopore simulators against seq2squiggle.

The benchmark includes four broad classes:
- Human read mode
- Human genome mode
- D.melanogaster genome mode
- Human variant calling

## Details

### Datasets

We use only R10.4.1 datasets from ONT for this benchmarking


|   Dataset	|   # of reads	|  Publication 	|   Accession ID	|   Description   	|
| -------- | ------- | -------- | ------- | ------- |
|D. melanogaster	| 229,767	| https://labs.epi2me.io/open-data-drosophila/  	| PRJNA914057 	| R10.4.1 Kit 14 (4kHz)	|
|   Human HCT116 	 |   5,149,102 	|   https://hasindu2008.github.io/f5c/docs/r10train 	|  PRJEB64592 	|   DNA Methylation standard. Methylated and non-methylated datasets are provided. Basecalled in study with Guppy 6.4.2 dna_r10.4.1_e8.2_400bps_sup	|


## Details 

### Training of seq2squiggle
For training, we utilize the non-methylated human dataset HCT116. 
The raw signal files are base-called using buttery-eel and aligned to the hg38 reference genome (excluding alternate contigs) with minimap2. 
Segmentation of the signal is performed using Uncalled4.

The dataset is split into training, validation, and test sets based on chromosomes. For the latest details, please refer to the configuration file.

### Read mode
In Read-mode, signals are simulated based on experimentally base-called reads (FASTQ/A data), allowing for direct comparison since both simulators use the same set of reads. 
This approach enables us to assess signal similarity and basecalling metrics. 

### Genome mode
In Genome-mode, the simulators generate reads based on the same input genome (FASTA data), resulting in different reads that should still resemble real reads in quality. This mode allows us to evaluate the overall basecalling accuracy of the simulators in replicating the characteristics of real nanopore sequencing reads.

### Variant calling
In variant calling mode, high-confidence NA12878 variants from Genome in a Bottle (v3.3.2) were integrated into the human reference genome. The simulators generate reads based on this modified diploid genome, which are then basecalled and aligned to the default hg38 reference genome. Variant calling is performed using Clair3 and the results were evaluated using RTGtools against the integrated high-confidence variants.


## Installation and running

The pipeline requires snakemake and mamba to run:

```
mamba env create -p ./env -f environment.yml
```

Next you can launch it on the cluster by submitting the slurm.sh submit script:

```
sbatch slurm.sh
```

Most output is written to the results/ directory, but benchmarking logs are written to the benchmarks/ directory.
You can customize the pipeline parameters by editing the config/configs.yml file or by modifying the rules defined in the workflow/rules/ directory.

### Using SLURM

The `workflow/slurm/` directory is based on Snakemake's SLURM profile, which allows our
snakemake job to automatically submit sub-jobs via SLURM. It was be installed
using cookiecutter:
```
    cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```
When asked if you want to set "SBATCH_DEFAULTS", the following was set:
```
    --output=logs/slurm/%x-%A_%a.out --error=logs/slurm/%x-%A_%a.err --open-mode=truncate
```
This can be seen in the file `slurm/settings.json`.
