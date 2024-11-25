# Do the actual benchmark runs of all software here


# Run simulators



rule run_squigulator_ideal_fly:
    input:
        config["d_melanogaster_ref"],
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
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule blow5_to_pod5_squigulator_ideal_fly:
    input:
        rules.run_squigulator_ideal_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_ideal.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """  

rule basecall_squigulator_ideal_fly:
    input:
        rules.blow5_to_pod5_squigulator_ideal_fly.output,
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
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
        ref=config["d_melanogaster_ref"]
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
        config["d_melanogaster_ref"],
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
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule blow5_to_pod5_squigulator_idealtime_fly:
    input:
        rules.run_squigulator_idealtime_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_idealtime.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_squigulator_idealtime_fly:
    input:
        rules.blow5_to_pod5_squigulator_idealtime_fly.output,
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
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
        ref=config["d_melanogaster_ref"]
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
        config["d_melanogaster_ref"],
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
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule blow5_to_pod5_squigulator_idealamp_fly:
    input:
        rules.run_squigulator_idealamp_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_idealamp.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_squigulator_idealamp_fly:
    input:
        rules.blow5_to_pod5_squigulator_idealamp_fly.output,
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
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
        ref=config["d_melanogaster_ref"]
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
        config["d_melanogaster_ref"],
    output:
        "results/d-melanogaster-noise/squigulator_reads_default.slow5",
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["squigulator_default"]
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
        resources/squigulator-v0.4.0/squigulator {params.commands} {input} -o {output} -n {params.sample_rate}
        """

rule blow5_to_pod5_squigulator_default_fly:
    input:
        rules.run_squigulator_default_fly.output,
    output:
        "results/d-melanogaster/squigulator_reads_default.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_squigulator_default_fly:
    input:
        rules.blow5_to_pod5_squigulator_default_fly.output,
    output:
        "results/d-melanogaster-noise/squigulator_reads_default.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
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
        ref=config["d_melanogaster_ref"]
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