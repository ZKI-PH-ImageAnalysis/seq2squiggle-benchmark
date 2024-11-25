# Seq2Squiggle Workflow with Wildcards for dwell and noise parameters
rule run_seq2squiggle_heatplot:
    input:
        model="results/train-human/seq2squiggle-logs/last.ckpt",
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise_heatmap"]["seq2squiggle"],
        dwell_std="{dwell_std}",
        noise_std="{noise_std}"
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_dwell-{dwell_std}_noise-{noise_std}.txt"
    log:
        "results/logs/run_seq2squiggle_dwell-{dwell_std}_noise-{noise_std}.log"
    threads: 32
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki,zki_mig"
    shell:
        """
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} \
            --model {input.model} -n {params.sample_rate} {params.commands} \
            --dwell-std {params.dwell_std} --noise-std {params.noise_std} \
            -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_heatplot:
    input:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads.pod5"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads.fastq"
    params:
        model=config["basecalling_model"]
    threads: 16
    log:
        "results/logs/basecall_seq2squiggle_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=30000,
        slurm="gpus=1",
        partition="zki"
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_heatplot:
    input:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads.fastq"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads_fixed.fastq"
    threads: 4
    log:
        "results/logs/remove_header_seq2squiggle_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=10000
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule align_seq2squiggle_heatplot:
    input:
        seqs="results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads_fixed.fastq",
        ref=config["d_melanogaster_ref"]
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/seq2squiggle_reads.sam"
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align_seq2squiggle_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=50000
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """