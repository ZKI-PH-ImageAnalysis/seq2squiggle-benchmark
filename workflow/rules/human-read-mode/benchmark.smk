# Do the actual benchmark runs of all software here


# Run simulators

rule run_seq2squiggle:
    input:
        model=rules.train_seq2squiggle.output,
        fastq=rules.remove_short_reads_subset.output,
        config="config/seq2squiggle-config.yml",
    output:
        "results/human/seq2squiggle_reads.slow5"
    params:
        commands=config["read_mode"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle.yml"
    log:
        "results/logs/seq2squiggle-simulation.log"
    threads: 16
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} {params.commands} -o {output} {input.fastq}
        """

rule fix_seq2squiggle_slow5:
    input:
        rules.run_seq2squiggle.output,
    output:
        "results/human/seq2squiggle_reads_fixed.slow5",
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

rule basecall_seq2squiggle:
    input:
        rules.run_seq2squiggle.output
    output:
        "results/human/seq2squiggle_reads.fastq"
    params:
        model=config["basecalling_model"]
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --use_tcp --port 8010 --slow5_threads 32
        """

rule remove_header_seq2squiggle_fastq:
    input:
        rules.basecall_seq2squiggle.output,
    output:
        "results/human/seq2squiggle_reads_fixed.fastq",
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

rule run_squigulator:
    input:
        rules.remove_short_reads_subset.output,
    output:
        "results/human/squigulator_reads.slow5",
    params:
        commands=config["read_mode"]["squigulator"],
    threads: 128
    log:
        "results/logs/squigulator-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        resources/squigulator/squigulator {params.commands} {input} -o {output}
        """

rule fix_squigulator_slow5:
    input:
        rules.run_squigulator.output,
    output:
        "results/human/squigulator_reads_fixed.slow5",
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

rule basecall_squigulator:
    input:
        rules.fix_squigulator_slow5.output,
    output:
        "results/human/squigulator_reads.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --use_tcp --port 8020 --slow5_threads 32
        """

rule remove_header_squigulator_fastq:
    input:
        rules.basecall_squigulator.output,
    output:
        "results/human/squigulator_reads_fixed.fastq",
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

rule align_squigulator:
    input:
        seqs=rules.remove_header_squigulator_fastq.output,
        ref="data/zymo-human/GRCh38_no_alt.fna"
    output:
        "results/human/squigulator_reads.sam",
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

rule align_seq2squiggle:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq.output,
        ref="data/zymo-human/GRCh38_no_alt.fna"
    output:
        "results/human/seq2squiggle_reads.sam",
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
