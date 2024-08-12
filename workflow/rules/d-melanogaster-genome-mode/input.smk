# Download genomes and pre-process by basecalling, assembly and nanopolish
import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

#localrules:
#    download_signal_data,
#    download_genome

#rule download_signal_data:
#    input:
#        #TODO,
#    output:
#        #TODO,
#    log:
#        "results/logs/runtime/download-data-d.log",
#    run:
#        shell("aws s3 sync --no-sign-request s3://ont-open-data/contrib/melanogaster_bkim_2023.01/flowcells/D.melanogaster.R1041.400bps/D_melanogaster_1/20221217_1251_MN20261_FAV70669_117da01a/ D_melanogaster_R1041_400bps_1")

#rule download_genome:
#    input:
#        #TODO
#    output:
#        #TODO
#    log:
#        "results/logs/runtime/download-human-genome.log",
#    run:
#        shell("mv {input} {output}")
#    # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz


rule basecall_fly_all:
    input:
        "data/d-melanogaster/file.slow5",
    output:
        "data/d-melanogaster/file.fastq",
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall-fly-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(3600),
        mem_mb=150000,
        slurm="gpus=2",
        partition="zki",   
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8020 --use_tcp --slow5_threads 16 --procs 32
        """

rule align_reference_all_fly:
    input:
        seqs=rules.basecall_fly_all.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "data/d-melanogaster/file.bam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/align-reference-fly-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """

rule sort_index_bam_fly:
    input:
        rules.align_reference_all_fly.output,
    output:
        "data/d-melanogaster/file_sorted.bam"
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_bam.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule test_bam_to_fastq_fly:
    input:
        rules.sort_index_bam_fly.output,
    output:
        "data/d-melanogaster/file_all.fastq"
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/bam_to_fastq.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools bam2fq {input} > {output} 
        """


rule subsample_slow5_data_fly:
    input:
        fastq=rules.test_bam_to_fastq_fly.output,
        slow5="data/d-melanogaster/file.slow5" 
    output:
        read_l="data/d-melanogaster/file_subsample.txt",
        slow5="data/d-melanogaster/file_subsample.slow5",
    params:
        sample_rate=config["subsample-fly"],
        header_file="data/d-melanogaster/header.slow5",
        tmp_slow5="data/d-melanogaster/PGXX22563_reads_tmp.slow5",
    threads: 128
    log:
        "results/logs/subsample-slow5-fly.log",
    resources:
        disk_mb=50000,
        mem_mb=100000,
    shell:
        # get read_list with x random headers from fastq and use slow5tools to extract these
        # removes @ symbol since this causes problems
        """
        awk '/^@/{{sub(/^@/, ""); print $1}}' {input.fastq} | shuf -n {params.sample_rate} > {output.read_l}
        resources/slow5tools-v1.1.0/slow5tools get -t {threads} {input.slow5} -l {output.read_l} -o {output.slow5}
        """


rule basecall_subset_fly:
    input:
        rules.subsample_slow5_data_fly.output.slow5,
    output:
        "data/d-melanogaster/file_subsample.fastq",
    params:
        model=config["basecalling_model"]
    threads: 128
    log:
        "results/logs/basecall_subset_fly.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,
        slurm="gpus=2",
        partition="zki",   
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8013 --use_tcp --slow5_threads 32 --procs 64
        """

rule remove_header_fastq_subset_fly:
    input:
        rules.basecall_subset_fly.output,
    output:
        "data/d-melanogaster/file_subsample_fixed.fastq",
    threads: 64
    log:
        "results/logs/fix-fastq-header.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule remove_short_reads_subset_fly:
    input:
        rules.remove_header_fastq_subset_fly.output,
    output:
        "data/d-melanogaster/file_subsample_fixed_2.fastq",
    params:
        max_len=config["max_seq_len"]
    threads: 64
    log:
        "results/logs/fix-fastq-header.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.max_len}) {{print header, seq, qheader, qseq}}}}' < {input} > {output}
        """

rule align_reference_subset_fly:
    input:
        seqs=rules.remove_short_reads_subset_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "data/d-melanogaster/file_subsample.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-reference.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """
