# Do the actual benchmark runs of all software here


# Run simulators

rule run_seq2squiggle_gm_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogasquiggle_reads.slow5"
    params:
        sample_rate=config["subsample-fly"],
        commands=config["d_melanogaster_genome_mode"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle.yml"
    benchmark:
        "results/benchmarks/seq2signal_fly.txt"
    log:
        "results/logs/run_seq2squiggle_gm_fly.log"
    threads: 32
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=150000,
        slurm="gpus=1",
        partition="zki,zki_mig", 
    shell:
        """
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_gm_fly:
    input:
        rules.run_seq2squiggle_gm_fly.output,
    output:
        "results/d-melanogaster/seq2squiggle_reads_fixed.slow5",
    threads: 8
    log:
        "results/logs/fix_seq2squiggle_slow5_gm_fly.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(50),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} /^[@#]/ {{print; next}} {{split($1, a, "!"); $1 = a[2]; print}}' {input} > {output}
        """

rule basecall_seq2squiggle_gm_fly:
    input:
        rules.run_seq2squiggle_gm_fly.output
    output:
        "results/d-melanogaster/seq2squiggle_reads.fastq"
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall_seq2squiggle_gm_fly.log"
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=2",
        partition="zki", 
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8070 --use_tcp --slow5_threads 32
        """

rule remove_header_seq2squiggle_fastq_gm_fly:
    input:
        rules.basecall_seq2squiggle_gm_fly.output,
    output:
        "results/d-melanogaster/seq2squiggle_reads_fixed.fastq",
    threads: 64
    log:
        "results/logs/remove_header_seq2squiggle_fastq_gm_fly.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule run_squigulator_gm_fly:
    input:
        "data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "results/d-melanogaster/squigulator_reads.slow5",
    params:
        sample_rate=config["subsample-fly"],
        commands=config["d_melanogaster_genome_mode"]["squigulator"]
    threads: 128
    benchmark:
        "results/benchmarks/squigulator_fly.txt"
    log:
        "results/logs/squigulator-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        resources/squigulator/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule fix_squigulator_slow5_gm_fly:
    input:
        rules.run_squigulator_gm_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_fixed.slow5",
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

rule basecall_squigulator_gm_fly:
    input:
        rules.fix_squigulator_slow5_gm_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads.fastq"
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall_squigulator_gm_fly.log",
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=100000,
        slurm="gpus=1",
        partition="zki,zki_mig",   
    shell:
        """
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8071 --use_tcp --slow5_threads 32
        """

rule remove_header_squigulator_fastq_gm_fly:
    input:
        rules.basecall_squigulator_gm_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_fixed.fastq",
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

rule align_squigulator_gm_fly:
    input:
        seqs=rules.remove_header_squigulator_fastq_gm_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster/squigulator_reads.sam",
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

rule align_seq2squiggle_gm_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_gm_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster/seq2squiggle_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align_seq2squiggle_gm_fly.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """
