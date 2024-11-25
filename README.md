# seq2squiggle-benchmark
Collection of scripts for benchmarking of seq2squiggle with other simulators and experimental data.

## Setup and running
This repo contains an automated workflow to benchmark nanopore simulators against seq2squiggle.

The benchmark includes four broad classes:
- Read mode
- Genome mode
- Noisy study
- Variant calling

## Details

### Datasets

We use publicly available R9.4.1 and R10.4.1 datasets from ONT for this benchmarking:

| Dataset                 | # of Total Reads | Chemistry                  | Website                                                    |
|-------------------------|------------------|----------------------------|------------------------------------------------------------|
| D. melanogaster         | 229,767          | R10.4.1 Kit 14 (4kHz)      | [Link](https://labs.epi2me.io/open-data-drosophila/)       |
| Human HCT116 (non-methylated) | 5,149,102    | R10.4.1 Kit 14 (4kHz)      | [Link](https://hasindu2008.github.io/f5c/docs/r10train)    |
| E. coli                 | 81,037           | R10.4.1 Kit 14 (5kHz)      | [Link](http://ftp.sra.ebi.ac.uk/vol1/run/ERR138/ERR13848445/) |
| Human NA12878           | 500,000          | R9.4.1 Kit 10 (4kHz)       | [Link](https://slow5.bioinf.science/na12878_prom_sub_slow5) |
| SARS-CoV-2 SP1          | 1,382,016        | R9.4.1 Kit 10 (4kHz)       | [Link](https://slow5.bioinf.science/SP1-raw-mapped)        |


## Details 

### Training of seq2squiggle
For training the R10.4.1 model, we utilize the non-methylated human dataset HCT116.
For training the R9.4.1 model, we utilize the human dataset NA12878.
The raw signal files are converted to POD5 format, base-called using dorado and aligned to the hg38 reference genome (excluding alternate contigs) with minimap2. 
Segmentation of the signal is performed using Uncalled4.

The dataset is split into training, validation, and test sets based on chromosomes. For the latest details, please refer to the configuration file.

### Read mode
In Read-mode, signals are simulated based on experimentally base-called reads (FASTQ/A data), allowing for direct comparison since both simulators use the same set of reads. This mode allows for direct comparisons between simulators, as both are applied to the same set of reads. By using this approach, it is possible to assess both signal similarity and basecalling performance metrics, ensuring a thorough evaluation of the simulation process.

### Genome mode
In Genome-mode, the simulators generate reads based on the same input genome (FASTA data), resulting in different reads that should still resemble real reads in quality. This mode allows us to evaluate the overall basecalling accuracy of the simulators in replicating the characteristics of real nanopore sequencing reads. 

### Noise study
The Noise Study compares seq2Squiggle with squigulator in Genome Mode under various noise conditions. To conduct this study, we generated a heatmap that illustrates the median match rate for combinations of different amplitude noise and dwell time noise levels. Seq2Squiggle was evaluated both with and without its noise sampler and duration sampler modules. When these modules are disabled, seq2Squiggle relies on static noise distributions that are similar to those used by squigulator. This comparison provides insight into how different noise parameters affect signal quality and basecalling performance.

### Variant calling
In variant calling mode, high-confidence NA12878 variants from Genome in a Bottle (v3.3.2) were integrated into the human reference genome. The simulators generate reads based on this modified diploid genome, which are then basecalled and aligned to the default hg38 reference genome. Variant calling is performed using Clair3 and the results were evaluated using RTGtools against the integrated high-confidence variants. This mode highlights the simulators’ ability to generate biologically meaningful signals that can be used for downstream variant analysis.

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
