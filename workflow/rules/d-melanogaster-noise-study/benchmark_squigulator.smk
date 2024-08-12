# Do the actual benchmark runs of all software here


# Run simulators



rule run_squigulator_ideal_fly:
    input:
        "data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "results/d-melanogaster-noise/squigulator_reads_ideal.slow5",
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["squigulator_ideal"]
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

rule fix_squigulator_slow5_ideal_fly:
    input:
        rules.run_squigulator_ideal_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_ideal.slow5",
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

rule basecall_squigulator_ideal_fly:
    input:
        rules.fix_squigulator_slow5_ideal_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_ideal.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8003 --use_tcp --slow5_threads 32
        """

rule remove_header_squigulator_fastq_ideal_fly:
    input:
        rules.basecall_squigulator_ideal_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_ideal.fastq",
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

rule align_squigulator_ideal_fly:
    input:
        seqs=rules.remove_header_squigulator_fastq_ideal_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster-noise/squigulator_reads_ideal.sam",
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



rule run_squigulator_idealtime_fly:
    input:
        "data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealtime.slow5",
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["squigulator_idealtime"]
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

rule fix_squigulator_slow5_idealtime_fly:
    input:
        rules.run_squigulator_idealtime_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_idealtime.slow5",
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

rule basecall_squigulator_idealtime_fly:
    input:
        rules.fix_squigulator_slow5_idealtime_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealtime.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8001 --use_tcp --slow5_threads 32
        """

rule remove_header_squigulator_fastq_idealtime_fly:
    input:
        rules.basecall_squigulator_idealtime_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_idealtime.fastq",
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

rule align_squigulator_idealtime_fly:
    input:
        seqs=rules.remove_header_squigulator_fastq_idealtime_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealtime.sam",
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



rule run_squigulator_idealamp_fly:
    input:
        "data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealamp.slow5",
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["squigulator_idealamp"]
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

rule fix_squigulator_slow5_idealamp_fly:
    input:
        rules.run_squigulator_idealamp_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_idealamp.slow5",
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

rule basecall_squigulator_idealamp_fly:
    input:
        rules.fix_squigulator_slow5_idealamp_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealamp.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config {params.model} --port 8007 --use_tcp --slow5_threads 32
        """

rule remove_header_squigulator_fastq_idealamp_fly:
    input:
        rules.basecall_squigulator_idealamp_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_idealamp.fastq",
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

rule align_squigulator_idealamp_fly:
    input:
        seqs=rules.remove_header_squigulator_fastq_idealamp_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster-noise/squigulator_reads_idealamp.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """


rule run_squigulator_default_fly:
    input:
        "data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
    output:
        "results/d-melanogaster-noise/squigulator_reads_default.slow5",
    params:
        sample_rate=config["subsample-fly-noise"],
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
        resources/squigulator/squigulator -x dna-r10-prom --sample-rate 4000.0 {input} -o {output} -n {params.sample_rate}
        """

rule fix_squigulator_slow5_default_fly:
    input:
        rules.run_squigulator_default_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_default.slow5",
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

rule basecall_squigulator_default_fly:
    input:
        rules.fix_squigulator_slow5_default_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_default.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8008 --use_tcp --slow5_threads 32
        """

rule remove_header_squigulator_fastq_default_fly:
    input:
        rules.basecall_squigulator_default_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_fixed_default.fastq",
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

rule align_squigulator_default_fly:
    input:
        seqs=rules.remove_header_squigulator_fastq_default_fly.output,
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    output:
        "results/d-melanogaster-noise/squigulator_reads_default.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """