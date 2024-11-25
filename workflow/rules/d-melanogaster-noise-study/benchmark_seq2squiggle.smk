rule run_seq2squiggle_nonoise_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_nonoise.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_nonoise"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_nonoise_fly:
    input:
        rules.run_seq2squiggle_nonoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_nonoise.fastq"
    params:
        model=config["basecalling_model"]
    threads: 64
    log:
        "results/logs/basecall_seq2squiggle_gm_fly.log"
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

rule remove_header_seq2squiggle_fastq_nonoise_fly:
    input:
        rules.basecall_seq2squiggle_nonoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_nonoise.fastq",
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

rule align_seq2squiggle_nonoise_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_nonoise_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_nonoise.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """



rule run_seq2squiggle_manualnoise_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manualnoise.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_manualnoise"],
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """


rule basecall_seq2squiggle_manualnoise_fly:
    input:
        rules.run_seq2squiggle_manualnoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manualnoise.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_manualnoise_fly:
    input:
        rules.basecall_seq2squiggle_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_manualnoise.fastq",
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

rule align_seq2squiggle_manualnoise_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_manualnoise_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manualnoise.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """



rule run_seq2squiggle_noisesampler_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_noisesampler.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_noisesampler"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly_nonoise.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_noisesampler_fly:
    input:
        rules.run_seq2squiggle_noisesampler_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_noisesampler.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_noisesampler_fly:
    input:
        rules.basecall_seq2squiggle_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_noisesampler.fastq",
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

rule align_seq2squiggle_noisesampler_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_noisesampler_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_noisesampler.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """

rule run_seq2squiggle_lengthregulator_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_lengthregulator_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_lengthregulator_fly:
    input:
        rules.basecall_seq2squiggle_lengthregulator_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator.fastq",
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

rule align_seq2squiggle_lengthregulator_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_lengthregulator_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """


rule run_seq2squiggle_lengthregulator_and_manualnoise_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_manualnoise.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator_manualnoise"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_lengthregulator_and_manualnoise_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_manualnoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_manualnoise.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_lengthregulator_and_manualnoise_fly:
    input:
        rules.basecall_seq2squiggle_lengthregulator_and_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator_and_manualnoise.fastq",
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

rule align_seq2squiggle_lengthregulator_and_manualnoise_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_lengthregulator_and_manualnoise_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_manualnoise.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """



rule run_seq2squiggle_manuallengthsampling_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manuallengthsampling.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_manuallength"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_manuallengthsampling_fly:
    input:
        rules.run_seq2squiggle_manuallengthsampling_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manuallengthsampling.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_manuallengthsampling_fly:
    input:
        rules.basecall_seq2squiggle_manuallengthsampling_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_manuallengthsampling.fastq",
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

rule align_seq2squiggle_manuallengthsampling_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_manuallengthsampling_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manuallengthsampling.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """






rule run_seq2squiggle_lengthregulator_and_noisesampler_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_noisesampler.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator_noisesampler"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule basecall_seq2squiggle_lengthregulator_and_noisesampler_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_noisesampler_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_noisesampler.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_lengthregulator_and_noisesampler_fly:
    input:
        rules.basecall_seq2squiggle_lengthregulator_and_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator_and_noisesampler.fastq",
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

rule align_seq2squiggle_lengthregulator_and_noisesampler_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_lengthregulator_and_noisesampler_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_noisesampler.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """


rule run_seq2squiggle_mLR_mNS_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta=config["d_melanogaster_ref"],
        config="config/seq2squiggle-config.yml",
        installed=rules.install_poetry.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_mLR_mNS.pod5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_mLR_mNS"]
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_fly.txt"
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
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """


rule basecall_seq2squiggle_mLR_mNS_fly:
    input:
        rules.run_seq2squiggle_mLR_mNS_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_mLR_mNS.fastq"
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
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_seq2squiggle_fastq_mLR_mNS_fly:
    input:
        rules.basecall_seq2squiggle_mLR_mNS_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_mLR_mNS.fastq",
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

rule align_seq2squiggle_mLR_mNS_fly:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_mLR_mNS_fly.output,
        ref=config["d_melanogaster_ref"]
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_mLR_mNS.sam",
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
        minimap2 -x map-ont --eqx -a --secondary=no -t {threads} {input.ref} {input.seqs} > {output}
        """
