# Do the actual benchmark runs of all software here


# Run simulators

rule run_seq2squiggle_gm_ecoli:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["e_coli_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/e-coli/seq2squiggle_reads.pod5"
    params:
        sample_rate=config["subsample-ecoli"],
        commands=config["e_coli_genome_mode"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_ecoli.txt"
    log:
        "results/logs/run_seq2squiggle_gm_ecoli.log"
    threads: 32
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_gm_ecoli:
    input:
        rules.run_seq2squiggle_gm_ecoli.output
    output:
        "results/e-coli/seq2squiggle_reads.fastq"
    params:
        model=config["basecalling_model_r10_5khz"]
    threads: 64
    log:
        "results/logs/basecall_seq2squiggle_gm_ecoli.log"
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki", 
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_gm_ecoli:
    input:
        rules.basecall_seq2squiggle_gm_ecoli.output,
    output:
        "results/e-coli/seq2squiggle_reads_fixed.fastq",
    threads: 64
    log:
        "results/logs/remove_header_seq2squiggle_fastq_gm_ecoli.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule run_squigulator_gm_ecoli:
    input:
        config["e_coli_ref"],
    output:
        "results/e-coli/squigulator_reads.blow5",
    params:
        sample_rate=config["subsample-ecoli"],
        commands=config["e_coli_genome_mode"]["squigulator"]
    threads: 128
    benchmark:
        "results/benchmarks/squigulator_ecoli.txt"
    log:
        "results/logs/squigulator-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule fix_squigulator_slow5_gm_ecoli:
    input:
        rules.run_squigulator_gm_ecoli.output,
    output:
        "results/e-coli/squigulator_reads_fixed.slow5",
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


rule blow5_to_pod5_squigulator_gm_ecoli:
    input:
        rules.run_squigulator_gm_ecoli.output,
    output:
        "results/e-coli/squigulator_reads_fixed.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """        


rule basecall_squigulator_gm_ecoli:
    input:
        rules.blow5_to_pod5_squigulator_gm_ecoli.output,
    output:
        "results/e-coli/squigulator_reads.fastq"
    params:
        model=config["basecalling_model_r10_5khz"]
    threads: 64
    log:
        "results/logs/basecall_squigulator_gm_ecoli.log",
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

rule remove_header_squigulator_fastq_gm_ecoli:
    input:
        rules.basecall_squigulator_gm_ecoli.output,
    output:
        "results/e-coli/squigulator_reads_fixed.fastq",
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

rule align_squigulator_gm_ecoli:
    input:
        seqs=rules.remove_header_squigulator_fastq_gm_ecoli.output,
        ref=config["e_coli_ref"]
    output:
        "results/e-coli/squigulator_reads.sam",
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

rule align_seq2squiggle_gm_ecoli:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_gm_ecoli.output,
        ref=config["e_coli_ref"]
    output:
        "results/e-coli/seq2squiggle_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align_seq2squiggle_gm_ecoli.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """
