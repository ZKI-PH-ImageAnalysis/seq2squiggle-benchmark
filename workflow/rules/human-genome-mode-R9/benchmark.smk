# Do the actual benchmark runs of all software here
# Run simulators

rule run_seq2squiggle_gm_R9:
    input:
        model=rules.train_seq2squiggle_R9.output,
        fasta=config["human_ref_chr22"],
        config="config/seq2squiggle-config-R9.yml",
        #installed=rules.install_poetry.output,
    output:
        "results/human_gm_R9/seq2squiggle_reads.pod5"
    params:
        sample_rate=config["subsample-human-R9"],
        seq2squiggle_params=config["genome_mode_R9"]["seq2squiggle"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_gm_R9.txt"
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

rule basecall_seq2squiggle_gm_R9:
    input:
        rules.run_seq2squiggle_gm_R9.output
    output:
        "results/human_gm_R9/seq2squiggle_reads.fastq"
    threads: 32
    params:
        model=config["basecalling_model_R9"]
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

rule remove_header_seq2squiggle_fastq_gm_R9:
    input:
        rules.basecall_seq2squiggle_gm_R9.output,
    output:
        "results/human_gm_R9/seq2squiggle_reads_fixed.fastq",
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

rule run_squigulator_gm_R9:
    input:
        config["human_ref_chr22"],
    output:
        "results/human_gm_R9/squigulator_reads.slow5",
    params:
        sample_rate=config["subsample-human-R9"],
        run_params=config["genome_mode_R9"]["squigulator"]
    threads: 64
    benchmark:
        "results/benchmarks/squigulator_gm_R9.txt"
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

rule fix_squigulator_slow5_gm_R9:
    input:
        rules.run_squigulator_gm_R9.output,
    output:
        "results/human_gm_R9/squigulator_reads_fixed.slow5",
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

rule blow5_to_pod5_squigulator_gm_R9:
    input:
        rules.run_squigulator_gm_R9.output,
    output:
        "results/human_gm_R9/squigulator_reads_fixed.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """        

rule basecall_squigulator_gm_R9:
    input:
        rules.blow5_to_pod5_squigulator_gm_R9.output,
    output:
        "results/human_gm_R9/squigulator_reads.fastq"
    params:
        model=config["basecalling_model_R9"]
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

rule remove_header_squigulator_fastq_gm_R9:
    input:
        rules.basecall_squigulator_gm_R9.output,
    output:
        "results/human_gm_R9/squigulator_reads_fixed.fastq",
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


rule run_deepsimulator_CI_gm_R9:
    input:
        config["human_ref_chr22"],
    output:
        directory("results/human_gm_R9/deepsimulator_CI_reads/raw_reads"),
    conda:
        "../../envs/deepsimulator.yml"
    params:
        sample_rate=config["subsample-human-R9"],
        run_params=config["genome_mode_R9"]["deepsimulator_CI"]
    threads: 64
    benchmark:
        "results/benchmarks/deepsimulator_CI_gm_R9.txt"
    log:
        "results/logs/deepsimulator_CI-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        ./resources/DeepSimulator_benchmark/deep_simulator.sh -i {input} -o {output} -n {params.sample_rate} {params.run_params} -c {threads}
        """


rule fast5_to_blow5_CI_gm_R9:
    input:
        rules.run_deepsimulator_CI_gm_R9.output,
    output:
        directory("results/human_gm_R9/deepsimulator_CI_reads/blow5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        ./resources/slow5tools-v1.1.0/slow5tools f2s {input}/fast5/ -d {output}
        """         

rule blow5_to_pod5_CI_gm_R9:
    input:
        rules.fast5_to_blow5_CI_gm_R9.output,
    output:
        directory("results/human_gm_R9/deepsimulator_CI_reads/pod5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -d {output}
        """  


rule basecall_deepsimulator_CI_gm_R9:
    input:
        rules.blow5_to_pod5_CI_gm_R9.output,
    output:
        "results/human_gm_R9/deepsimulator_CI_reads.fastq"
    params:
        model=config["basecalling_model_R9"]
    threads: 64
    log:
        "results/logs/basecall-human-deepsimulator_CI.log",
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



rule align_deepsimulator_CI_gm_R9:
    input:
        seqs=rules.basecall_deepsimulator_CI_gm_R9.output,
        ref=config["human_ref"],
    output:
        "results/human_gm_R9/deepsimulator_CI_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-deepsimulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """



rule run_deepsimulator_CD_gm_R9:
    input:
        config["human_ref_chr22"],
    output:
        directory("results/human_gm_R9/deepsimulator_CD_reads/raw_reads"),
    conda:
        "../../envs/deepsimulator.yml"
    params:
        sample_rate=config["subsample-human-R9"],
        run_params=config["genome_mode_R9"]["deepsimulator_CD"]
    threads: 64
    benchmark:
        "results/benchmarks/deepsimulator_CD_gm_R9.txt"
    log:
        "results/logs/deepsimulator_CD-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        ./resources/DeepSimulator_benchmark/deep_simulator.sh -i {input} -o {output} -n {params.sample_rate} {params.run_params} -c {threads}
        """


rule fast5_to_blow5_CD_gm_R9:
    input:
        rules.run_deepsimulator_CD_gm_R9.output,
    output:
        directory("results/human_gm_R9/deepsimulator_CD_reads/blow5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        ./resources/slow5tools-v1.1.0/slow5tools f2s {input}/fast5/ -d {output}
        """         

rule blow5_to_pod5_CD_gm_R9:
    input:
        rules.fast5_to_blow5_CD_gm_R9.output,
    output:
        directory("results/human_gm_R9/deepsimulator_CD_reads/pod5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -d {output}
        """  


rule basecall_deepsimulator_CD_gm_R9:
    input:
        rules.blow5_to_pod5_CD_gm_R9.output,
    output:
        "results/human_gm_R9/deepsimulator_CD_reads.fastq"
    params:
        model=config["basecalling_model_R9"]
    threads: 64
    log:
        "results/logs/basecall-human-deepsimulator_CD.log",
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



rule align_deepsimulator_CD_gm_R9:
    input:
        seqs=rules.basecall_deepsimulator_CD_gm_R9.output,
        ref=config["human_ref"],
    output:
        "results/human_gm_R9/deepsimulator_CD_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-deepsimulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """





rule align_squigulator_gm_R9:
    input:
        seqs=rules.remove_header_squigulator_fastq_gm_R9.output,
        ref=config["human_ref"],
    output:
        "results/human_gm_R9/squigulator_reads.sam",
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

rule align_seq2squiggle_gm_R9:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_gm_R9.output,
        ref=config["human_ref"],
    output:
        "results/human_gm_R9/seq2squiggle_reads.sam",
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
