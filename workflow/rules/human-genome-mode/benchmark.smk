# Do the actual benchmark runs of all software here
# Run simulators



rule run_seq2squiggle_gm:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["human_ref_chr22"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/human_gm/seq2squiggle_reads.pod5"
    params:
        sample_rate=config["subsample-human"],
        seq2squiggle_params=config["genome_mode"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_gm.txt"
    threads: 64
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """ 
        poetry -C resources/seq2squiggle/ run seq2squiggle predict {params.seq2squiggle_params} --config {input.config} --model {input.model} -n {params.sample_rate} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_gm:
    input:
        rules.run_seq2squiggle_gm.output
    output:
        "results/human_gm/seq2squiggle_reads.fastq"
    threads: 32
    params:
        model=config["basecalling_model"]
    log:
        "results/logs/basecall-human-seq2squiggle.log"
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_gm:
    input:
        rules.basecall_seq2squiggle_gm.output,
    output:
        "results/human_gm/seq2squiggle_reads_fixed.fastq",
    threads: 64
    log:
        "results/logs/fix-reads-seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule run_squigulator_gm:
    input:
        config["human_ref_chr22"],
    output:
        "results/human_gm/squigulator_reads.slow5",
    params:
        sample_rate=config["subsample-human"],
        run_params=config["genome_mode"]["squigulator"]
    threads: 64
    benchmark:
        "results/benchmarks/squigulator.txt"
    log:
        "results/logs/squigulator-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        resources/squigulator-v0.4.0/squigulator {params.run_params} {input} -o {output} -n {params.sample_rate}
        """

rule fix_squigulator_slow5_gm:
    input:
        rules.fix_squigulator_slow5_rm_readids.output,
    output:
        "results/human_gm/squigulator_reads_fixed.slow5",
    threads: 8
    log:
        "results/logs/squigulator-fixslow5.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(50),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} /^[@#]/ {{print; next}} {{split($1, a, "!"); $1 = a[2]; print}}' {input} > {output}
        """


rule blow5_to_pod5_squigulator_gm:
    input:
        rules.run_squigulator_gm.output,
    output:
        "results/human_gm/squigulator_reads_fixed.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """        

rule basecall_squigulator_gm:
    input:
        rules.blow5_to_pod5_squigulator_gm.output,
    output:
        "results/human_gm/squigulator_reads.fastq"
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall-human-squigulator.log",
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki,zki_mig",   
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_squigulator_fastq_gm:
    input:
        rules.basecall_squigulator_gm.output,
    output:
        "results/human_gm/squigulator_reads_fixed.fastq",
    threads: 64
    log:
        "results/logs/fix-reads-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule align_squigulator_gm:
    input:
        seqs=rules.remove_header_squigulator_fastq_gm.output,
        ref=config["human_ref"],
    output:
        "results/human_gm/squigulator_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """

rule align_seq2squiggle_gm:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_gm.output,
        ref=config["human_ref"],
    output:
        "results/human_gm/seq2squiggle_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """
