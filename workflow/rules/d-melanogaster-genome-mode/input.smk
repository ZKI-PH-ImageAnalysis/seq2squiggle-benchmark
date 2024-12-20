# Download genomes and pre-process by basecalling, assembly and nanopolish

rule blow5_to_pod5_fly:
    input:
        config["d_melanogaster_signal"],
    output:
        "data/processed-data/d-melanogaster/reads.pod5"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """


rule basecall_fly_all:
    input:
        rules.blow5_to_pod5_fly.output,
    output:
        "data/processed-data/d-melanogaster/reads.fastq",
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall-fly-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(3600),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki",   
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule align_reference_all_fly:
    input:
        seqs=rules.basecall_fly_all.output,
        ref=config["d_melanogaster_ref"],
    output:
        "data/processed-data/d-melanogaster/reads.bam",
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
        "data/processed-data/d-melanogaster/reads_sorted.bam"
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
        "data/processed-data/d-melanogaster/file_all.fastq"
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
        slow5=config["d_melanogaster_signal"] 
    output:
        read_l="data/processed-data/d-melanogaster/reads_subset.txt",
        slow5="data/processed-data/d-melanogaster/reads_subset.slow5",
    params:
        sample_rate=config["subsample-fly"],
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
        awk 'NR % 4 == 1 {{sub(/^@/, ""); print $1}}' {input.fastq} | shuf -n {params.sample_rate} > {output.read_l}
        resources/slow5tools-v1.1.0/slow5tools get -t {threads} {input.slow5} -l {output.read_l} -o {output.slow5} --skip
        resources/slow5tools-v1.1.0/slow5tools index {output.slow5}
        """


rule blow5_to_pod5_fly_subset:
    input:
        rules.subsample_slow5_data_fly.output.slow5,
    output:
        "data/processed-data/d-melanogaster/reads_subset.pod5"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """


rule basecall_subset_fly:
    input:
        rules.blow5_to_pod5_fly_subset.output,
    output:
        "data/processed-data/d-melanogaster/reads_subsample_pass.fastq",
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_fastq_subset_fly:
    input:
        rules.basecall_subset_fly.output,
    output:
        "data/processed-data/d-melanogaster/file_subsample_fixed.fastq",
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
        "data/processed-data/d-melanogaster/file_subsample_fixed_2.fastq",
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
        ref=config["d_melanogaster_ref"],
    output:
        "data/processed-data/d-melanogaster/file_subsample_pass.sam",
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
