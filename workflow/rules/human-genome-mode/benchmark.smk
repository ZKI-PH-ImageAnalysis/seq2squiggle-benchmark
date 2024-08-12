# Do the actual benchmark runs of all software here


# Run simulators

rule run_seq2squiggle_gm:
    input:
        model=rules.train_seq2squiggle.output,
        fasta="data/zymo-human/GRCh38_no_alt.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/human_gm/seq2squiggle_reads.slow5"
    params:
        sample_rate=config["subsample-human"],
        seq2squiggle_params=config["genome_mode"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle.yml"
    benchmark:
        "results/benchmarks/seq2squiggle.txt"
    log:
        "results/logs/seq2squiggle-simulation.log"
    threads: 64
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict {params.seq2squiggle_params} --config {input.config} --model {input.model} -n {params.sample_rate} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_gm:
    input:
        rules.run_seq2squiggle_gm.output,
    output:
        "results/human_gm/seq2squiggle_reads_fixed.slow5",
    threads: 8
    log:
        "results/logs/seq2squiggle-fixslow5.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(50),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} /^[@#]/ {{print; next}} {{split($1, a, "!"); $1 = a[2]; print}}' {input} > {output}
        """

rule basecall_seq2squiggle_gm:
    input:
        rules.run_seq2squiggle_gm.output
    output:
        "results/human_gm/seq2squiggle_reads.fastq"
    threads: 32
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --use_tcp --port 8017 --slow5_threads 32
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
        "data/zymo-human/GRCh38_no_alt.fna",
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
        resources/squigulator/squigulator {params.run_params} {input} -o {output} -n {params.sample_rate}
        """

rule fix_squigulator_slow5_gm:
    input:
        rules.run_squigulator_gm.output,
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

rule basecall_squigulator_gm:
    input:
        rules.fix_squigulator_slow5_gm.output,
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --use_tcp --port 8016 --slow5_threads 32
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
        ref="data/zymo-human/GRCh38_no_alt.fna"
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
        ref="data/zymo-human/GRCh38_no_alt.fna"
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
