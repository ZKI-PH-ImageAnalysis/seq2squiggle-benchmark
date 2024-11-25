# Download genomes and pre-process by basecalling, assembly and nanopolish

rule blow5_to_pod5_sarscov_R9:
    input:
        config["sars_cov_signal"],
    output:
        "data/processed-data/sars-cov-R9/reads.pod5"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -o {output}
        """


rule basecall_sarscov_R9_all:
    input:
        rules.blow5_to_pod5_sarscov_R9.output,
    output:
        "data/processed-data/sars-cov-R9/reads.fastq",
    params:
        model=config["basecalling_model_R9"]
    threads: 64
    log:
        "results/logs/basecall-sarscov_R9-all.log",
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

rule align_reference_all_sarscov_R9:
    input:
        seqs=rules.basecall_sarscov_R9_all.output,
        ref=config["sars_cov_ref"],
    output:
        "data/processed-data/sars-cov-R9/reads.bam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/align-reference-sarscov_R9-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """

rule sort_index_bam_sarscov_R9:
    input:
        rules.align_reference_all_sarscov_R9.output,
    output:
        "data/processed-data/sars-cov-R9/reads_sorted.bam"
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

rule test_bam_to_fastq_sarscov_R9:
    input:
        rules.sort_index_bam_sarscov_R9.output,
    output:
        "data/processed-data/sars-cov-R9/file_all.fastq"
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


rule subsample_slow5_data_sarscov_R9:
    input:
        fastq=rules.test_bam_to_fastq_sarscov_R9.output,
        slow5=config["sars_cov_signal"] 
    output:
        read_l="data/processed-data/sars-cov-R9/reads_subset.txt",
        slow5="data/processed-data/sars-cov-R9/reads_subset.slow5",
    params:
        sample_rate=config["subsample-sarscov_R9"],
    threads: 128
    log:
        "results/logs/subsample-slow5-sarscov_R9.log",
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


rule blow5_to_pod5_sarscov_R9_subset:
    input:
        rules.subsample_slow5_data_sarscov_R9.output.slow5,
    output:
        "data/processed-data/sars-cov-R9/reads_subset.pod5"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """


rule basecall_subset_sarscov_R9:
    input:
        rules.blow5_to_pod5_sarscov_R9_subset.output,
    output:
        "data/processed-data/sars-cov-R9/reads_subsample_pass.fastq",
    params:
        model=config["basecalling_model_R9"]
    threads: 128
    log:
        "results/logs/basecall_subset_sarscov_R9.log",
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

rule remove_header_fastq_subset_sarscov_R9:
    input:
        rules.basecall_subset_sarscov_R9.output,
    output:
        "data/processed-data/sars-cov-R9/file_subsample_fixed.fastq",
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

rule remove_short_reads_subset_sarscov_R9:
    input:
        rules.remove_header_fastq_subset_sarscov_R9.output,
    output:
        "data/processed-data/sars-cov-R9/file_subsample_fixed_2.fastq",
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

rule align_reference_subset_sarscov_R9:
    input:
        seqs=rules.remove_short_reads_subset_sarscov_R9.output,
        ref=config["sars_cov_ref"],
    output:
        "data/processed-data/sars-cov-R9/file_subsample_pass.sam",
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
