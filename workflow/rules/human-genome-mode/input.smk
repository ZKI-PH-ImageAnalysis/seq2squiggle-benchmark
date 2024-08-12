# Download genomes and pre-process by basecalling, assembly and nanopolish
import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

#localrules:
#    download_signal_data,
#    download_genome

#rule download_signal_data:
#    input:
#        HTTP.remote(config["human_reference_genome"]),
#    output:
#        "data/zymo-human/PGXX22563_pcr/PGXX22563_reads.slow5",
#    log:
#        "results/logs/runtime/download-data.log",
#    run:
#        shell("mv {input} {output}")

#rule download_genome:
#    input:
#        HTTP.remote(config["human_signal_data"]),
#    output:
#        "data/zymo-human/GRCh38_no_alt.fna",
#    log:
#        "results/logs/runtime/download-human-genome.log",
#    run:
#        shell("mv {input} {output}")


rule basecall_human_all_uncalled4:
    input:
        "data/zymo-human/PGXX22563_pcr/PGXX22563_reads.slow5",
    output:
        "data/zymo-human/PGXX22563_all.sam",
    params:
        model=config["basecalling_model"]
    threads: 128
    log:
        "results/logs/basecall-human-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=100000,
        slurm="gpus=4",
        partition="zki",   
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8019 --use_tcp --slow5_threads 32 --procs 64 --moves_out
        """


rule sam_to_fastq_all:
    input:
        rules.basecall_human_all_uncalled4.output
    output:
        "data/zymo-human/PGXX22563_all.fastq",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sam_to_fastq_all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=200000,
    shell:
        """
        samtools fastq -T "mv,ts,pi,sp,ns" {input} > {output} 
        """


rule align_reference_all_uncalled4:
    input:
        seqs=rules.sam_to_fastq_all.output,
        ref="data/zymo-human/GRCh38_no_alt.fna",
    output:
        "data/zymo-human/PGXX22563_all.bam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/align-reference-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=250000,   
    shell:
        """
        minimap2 -y {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """

rule sort_index_all_uncalled4:
    input:
        rules.align_reference_all_uncalled4.output,
    output:
        "data/zymo-human/PGXX22563_all_sorted.bam"
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_bam.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=250000,   
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """


rule samtools_split_chr_uncalled4:
    input:
        rules.sort_index_all_uncalled4.output,
    output:
        train="data/zymo-human/PGXX22563_reads_train.bam",
        valid="data/zymo-human/PGXX22563_reads_valid.bam",
        test="data/zymo-human/PGXX22563_reads_test.bam",
    params:
        region_train=config["region_train"],
        region_valid=config["region_valid"],
        region_test=config["region_test"],
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/samtools_split_chr.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools view -@ {threads} {input} {params.region_train} -h -o {output.train}
        samtools view -@ {threads} {input} {params.region_valid} -h -o {output.valid}
        samtools view -@ {threads} {input} {params.region_test} -h -o {output.test}
        """


rule sort_index_train:
    input:
        rules.samtools_split_chr_uncalled4.output.train,
    output:
        "data/zymo-human/PGXX22563_reads_train_sorted.bam"
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


rule sort_index_valid:
    input:
        rules.samtools_split_chr_uncalled4.output.valid,
    output:
        "data/zymo-human/PGXX22563_reads_valid_sorted.bam"
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


rule uncalled4_eventalign_train:
    input:
        bam=rules.sort_index_train.output,
        ref="data/zymo-human/GRCh38_no_alt.fna",
        slow5="data/zymo-human/PGXX22563_pcr/PGXX22563_reads.slow5",
    output:
        "data/zymo-human/events-train.tsv"
    params:
        commands=config["uncalled4_commands"]
    threads: 128
    log: 
        "results/logs/uncalled4-align.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=100000,
        partition="zki",   
    shell:
        """
        uncalled4 align {input.ref} {input.slow5} --bam-in {input.bam} --eventalign-out {output} {params.commands}
        """


rule uncalled4_eventalign_valid:
    input:
        bam=rules.sort_index_valid.output,
        ref="data/zymo-human/GRCh38_no_alt.fna",
        slow5="data/zymo-human/PGXX22563_pcr/PGXX22563_reads.slow5",
    output:
        "data/zymo-human/events-valid.tsv"
    params:
        commands=config["uncalled4_commands"]
    threads: 128
    log: 
        "results/logs/uncalled4-align.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=100000,
        partition="zki",   
    shell:
        """
        uncalled4 align {input.ref} {input.slow5} --bam-in {input.bam} --eventalign-out {output} {params.commands}
        """


rule eventalign_train_norm:
    input:
        rules.uncalled4_eventalign_train.output
    output:
        "data/zymo-human/events-train-norm.tsv"
    conda:
        "../../envs/seq2squiggle.yml"
    threads: 128
    log: 
        "results/logs/uncalled4-norm.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=500000,
        partition="zki", 
    shell:
        """
        python workflow/scripts/standardize-events.py {input} {output}
        """


rule eventalign_valid_norm:
    input:
        rules.uncalled4_eventalign_valid.output
    output:
        "data/zymo-human/events-valid-norm.tsv"
    conda:
        "../../envs/seq2squiggle.yml"
    threads: 128
    log: 
        "results/logs/uncalled4-norm-valid.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=500000,
        partition="zki", 
    shell:
        """
        python workflow/scripts/standardize-events.py {input} {output}
        """

rule preprocess_train_seq2squiggle:
    input:
        train=rules.eventalign_train_norm.output,
        config="config/seq2squiggle-config.yml",
    output:
        directory("results/zymo-human/seq2squiggle-train"),
    conda:
        "../../envs/seq2squiggle.yml"
    threads: 128
    log:
        "results/logs/preprocess_seq2squiggle.log",
    resources:
        disk_mb=100000,
        runtime=add_slack(3600),
        mem_mb=350000,
        partition="zki"
    shell:
        """
        python resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py preprocess {input.train} {output} --config {input.config}
        """
        # --max-chunks 100000000 are ~ 3mil chunks if mapping is 32 to 500


rule preprocess_valid_seq2squiggle:
    input:
        valid=rules.eventalign_valid_norm.output,
        config="config/seq2squiggle-config.yml",
    output:
        directory("results/zymo-human/seq2squiggle-valid"),
    conda:
        "../../envs/seq2squiggle.yml"
    threads: 128
    log:
        "results/logs/preprocess_seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(3600),
        mem_mb=100000,
        partition="zki"
    shell:
        """
        python resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py preprocess {input.valid} {output} --config {input.config}
        """


rule train_seq2squiggle:
    input:
        train_dir=rules.preprocess_train_seq2squiggle.output,
        valid_dir=rules.preprocess_valid_seq2squiggle.output,
        config="config/seq2squiggle-config.yml",
    output:
        "results/zymo-human/seq2squiggle-logs/last.ckpt",
    conda:
        "../../envs/seq2squiggle.yml"
    threads: 64
    log:
        "results/logs/train_seq2squiggle.log",
    resources:
        disk_mb=50000, 
        runtime=add_slack(9000),
        mem_mb=150000, #300000
        slurm="gpus=8",
        partition="zki",     
    shell:
        """
        export HTTPS_PROXY="http://fw-bln.rki.local:8020"
        export https_proxy="http://fw-bln.rki.local:8020"
        python resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py train {input.train_dir} {input.valid_dir} \
         --config {input.config} --model {output}
        """


rule test_bam_to_fastq:
    input:
        rules.samtools_split_chr_uncalled4.output.test,
    output:
        "data/zymo-human/PGXX22563_reads_test.fastq"
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


rule subsample_slow5_data:
    input:
        fastq=rules.test_bam_to_fastq.output,
        slow5="data/zymo-human/PGXX22563_pcr/PGXX22563_reads.slow5" 
    output:
        read_l="data/zymo-human/PGXX22563_pcr/PGXX22563_reads_subsample.txt",
        slow5="data/zymo-human/PGXX22563_pcr/PGXX22563_reads_subsample.slow5",
    params:
        sample_rate=config["subsample-human"],
        header_file="data/zymo-human/PGXX22563_pcr/header.slow5",
        tmp_slow5="data/zymo-human/PGXX22563_pcr/PGXX22563_reads_tmp.slow5",
    threads: 128
    log:
        "results/logs/subsample-slow5-human.log",
    resources:
        disk_mb=50000,
        mem_mb=100000,
    shell:
        # get read_list with x random headers from fastq and use slow5tools to extract these
        # removes @ symbol since this causes problems
        """
        awk '/^@/{{sub(/^@/, ""); print $1}}' {input.fastq} | shuf -n {params.sample_rate} > {output.read_l}
        resources/slow5tools-v1.1.0/slow5tools get -t {threads} {input.slow5} -l {output.read_l} -o {output.slow5}
        resources/slow5tools-v1.1.0/slow5tools index {output.slow5}
        """


rule basecall_human_subset:
    input:
        rules.subsample_slow5_data.output.slow5,
    output:
        "data/zymo-human/reads_subsample_pass.fastq",
    params:
        model=config["basecalling_model"]
    threads: 128
    log:
        "results/logs/basecall_human_subset.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki",   
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --use_tcp --port 8020 --slow5_threads 32 --procs 64
        """

rule remove_header_fastq_subset:
    input:
        rules.basecall_human_subset.output,
    output:
        "data/zymo-human/reads_subsample_pass_fixed.fastq",
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

rule remove_short_reads_subset:
    input:
        rules.remove_header_fastq_subset.output,
    output:
        "data/zymo-human/reads_subsample_pass_fixed_2.fastq",
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

rule align_reference_subset:
    input:
        seqs=rules.remove_short_reads_subset.output,
        ref="data/zymo-human/GRCh38_no_alt.fna",
    output:
        "data/zymo-human/reads_subsample_pass.sam",
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
