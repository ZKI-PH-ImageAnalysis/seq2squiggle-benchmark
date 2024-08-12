rule run_seq2squiggle_nonoise_fly:
    input:
        model=rules.train_seq2squiggle.output,
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_nonoise.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_nonoise"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """


rule fix_seq2squiggle_slow5_nonoise_fly:
    input:
        rules.run_seq2squiggle_nonoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_nonoise.slow5",
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

rule basecall_seq2squiggle_nonoise_fly:
    input:
        rules.run_seq2squiggle_nonoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_nonoise.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8002 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manualnoise.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_manualnoise"],
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_manualnoise_fly:
    input:
        rules.run_seq2squiggle_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_manualnoise.slow5",
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

rule basecall_seq2squiggle_manualnoise_fly:
    input:
        rules.run_seq2squiggle_manualnoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manualnoise.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8043 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_noisesampler.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_noisesampler"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_noisesampler_fly:
    input:
        rules.run_seq2squiggle_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_noisesampler.slow5",
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

rule basecall_seq2squiggle_noisesampler_fly:
    input:
        rules.run_seq2squiggle_noisesampler_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_noisesampler.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8009 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """


rule fix_seq2squiggle_slow5_lengthregulator_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator.slow5",
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

rule basecall_seq2squiggle_lengthregulator_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8022 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_manualnoise.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator_manualnoise"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_lengthregulator_and_manualnoise_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator_and_manualnoise.slow5",
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

rule basecall_seq2squiggle_lengthregulator_and_manualnoise_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_manualnoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_manualnoise.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8025 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manuallengthsampling.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_manuallength"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """


rule fix_seq2squiggle_slow5_manuallengthsampling_fly:
    input:
        rules.run_seq2squiggle_manuallengthsampling_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_manuallengthsampling.slow5",
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

rule basecall_seq2squiggle_manuallengthsampling_fly:
    input:
        rules.run_seq2squiggle_manuallengthsampling_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_manuallengthsampling.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8024 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_noisesampler.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_lengthregulator_noisesampler"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_lengthregulator_and_noisesampler_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_lengthregulator_and_noisesampler.slow5",
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

rule basecall_seq2squiggle_lengthregulator_and_noisesampler_fly:
    input:
        rules.run_seq2squiggle_lengthregulator_and_noisesampler_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_lengthregulator_and_noisesampler.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8081 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
        fasta="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        config="config/seq2squiggle-config.yml",
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_mLR_mNS.slow5"
    params:
        sample_rate=config["subsample-fly-noise"],
        commands=config["d_melanogaster_noise"]["seq2squiggle_mLR_mNS"]
    conda:
        "../../envs/seq2squiggle.yml"
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
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --config {input.config} --model {input.model} -n {params.sample_rate} {params.commands} -o {output} {input.fasta}
        """

rule fix_seq2squiggle_slow5_mLR_mNS_fly:
    input:
        rules.run_seq2squiggle_mLR_mNS_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_fixed_mLR_mNS.slow5",
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

rule basecall_seq2squiggle_mLR_mNS_fly:
    input:
        rules.run_seq2squiggle_mLR_mNS_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle_reads_mLR_mNS.fastq"
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
        buttery-eel -i {input} -o {output} -x \"cuda:all\" -g resources/ont-guppy-6_5_7/bin/ --config dna_r10.4.1_e8.2_400bps_sup.cfg --port 8082 --use_tcp --slow5_threads 32
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
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
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
