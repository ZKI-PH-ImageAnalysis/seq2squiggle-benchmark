# Squigulator Workflow with Wildcards for dwell and noise parameters
rule run_squigulator_heatplot:
    input:
        config["d_melanogaster_ref"]
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise_heatmap"]["squigulator"],
        dwell_std="{dwell_std}",
        noise_std="{noise_std}"
    threads: 128
    benchmark:
        "results/benchmarks/squigulator_dwell-{dwell_std}_noise-{noise_std}.txt"
    log:
        "results/logs/squigulator-simulation_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000
    shell:
        """
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} \
            --dwell-std {params.dwell_std} --amp-noise {params.noise_std} \
            -o {output} -n {params.sample_rate}
        """

rule blow5_to_pod5_squigulator_heatplot:
    input:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.slow5"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_squigulator_heatplot:
    input:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.pod5"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.fastq"
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall_squigulator_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki,zki_mig"
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_squigulator_fastq_heatplot:
    input:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.fastq"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads_fixed.fastq"
    threads: 64
    log:
        "results/logs/fix-reads-squigulator_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule align_squigulator_heatplot:
    input:
        seqs="results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads_fixed.fastq",
        ref=config["d_melanogaster_ref"]
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/squigulator_reads.sam"
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-squigulator_dwell-{dwell_std}_noise-{noise_std}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """