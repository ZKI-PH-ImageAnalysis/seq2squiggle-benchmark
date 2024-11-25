# Download
# ./slow5curl reads https://gtgseq.s3.amazonaws.com/ont-r10-5khz-dna/NA24385/raw/PGXXXX230339_reads.blow5 > allreads.txt
# ./slow5curl get --listallreads.txt https://gtgseq.s3.amazonaws.com/ont-r10-5khz-dna/NA24385/raw/PGXXXX230339_reads.blow5 -o selected.blow5

rule blow5_to_pod5_human_all_R9:
    input:
        config["human_signal_R9"],
    output:
        "data/processed-data/R9/reads.pod5"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_human_all_uncalled4_R9:
    input:
        rules.blow5_to_pod5_human_all_R9.output,
    output:
        "data/processed-data/R9/reads.sam",
    params:
        model=config["basecalling_model_R9"]
    threads: 64
    log:
        "results/logs/basecall-human-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=50000,
        slurm="gpus=1",
        partition="zki,main",   
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-moves {params.model} {input} > {output}
        """

# TODO Filter reads below Q-Score of 9?

rule sam_to_fastq_all_R9:
    input:
        rules.basecall_human_all_uncalled4_R9.output
    output:
        "data/processed-data/R9/reads.fastq",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sam_to_fastq_all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=200000,
    shell:
        """
        samtools fastq -T "mv,ts,pi,sp,ns" {input} > {output} 
        """


rule align_reference_all_uncalled4_R9:
    input:
        seqs=rules.sam_to_fastq_all_R9.output,
        ref=config["human_ref"],
    output:
        "data/processed-data/R9/reads.bam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/align-reference-all.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=250000,   
    shell:
        """
        minimap2 -y {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """

rule sort_index_all_uncalled4_R9:
    input:
        rules.align_reference_all_uncalled4_R9.output,
    output:
        "data/processed-data/R9/reads_all_sorted.bam"
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_bam.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=250000,   
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """


rule samtools_split_chr_uncalled4_R9:
    input:
        rules.sort_index_all_uncalled4_R9.output,
    output:
        train="data/processed-data/R9/reads_train.bam",
        valid="data/processed-data/R9/reads_valid.bam",
        test="data/processed-data/R9/reads_test.bam",
    params:
        region_train=config["region_train"],
        region_valid=config["region_valid"],
        region_test=config["region_test"],
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/samtools_split_chr.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools view -@ {threads} {input} {params.region_train} -h -o {output.train}
        samtools view -@ {threads} {input} {params.region_valid} -h -o {output.valid}
        samtools view -@ {threads} {input} {params.region_test} -h -o {output.test}
        """


rule sort_index_train_R9:
    input:
        rules.samtools_split_chr_uncalled4_R9.output.train,
    output:
        "data/processed-data/R9/reads_train_sorted.bam"
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_bam.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule sort_index_valid_R9:
    input:
        rules.samtools_split_chr_uncalled4_R9.output.valid,
    output:
        "data/processed-data/R9/reads_valid_sorted.bam"
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_bam.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """


rule uncalled4_eventalign_train_R9:
    input:
        bam=rules.sort_index_train_R9.output,
        ref=config["human_ref"],
        slow5=config["human_signal_R9"],
    output:
        "data/processed-data/R9/events-train.tsv"
    params:
        commands=config["uncalled4_commands_R9"]
    threads: 128
    log: 
        "results/logs/uncalled4-align.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=100000,
        partition="zki",   
    shell:
        """
        uncalled4 align {input.ref} {input.slow5} --bam-in {input.bam} --eventalign-out {output} {params.commands}
        """


rule uncalled4_eventalign_valid_R9:
    input:
        bam=rules.sort_index_valid_R9.output,
        ref=config["human_ref"],
        slow5=config["human_signal_R9"],
    output:
        "data/processed-data/R9/events-valid.tsv"
    params:
        commands=config["uncalled4_commands_R9"]
    threads: 128
    log: 
        "results/logs/uncalled4-align.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=100000,
        partition="zki",   
    shell:
        """
        uncalled4 align {input.ref} {input.slow5} --bam-in {input.bam} --eventalign-out {output} {params.commands}
        """


rule eventalign_train_norm_R9:
    input:
        rules.uncalled4_eventalign_train_R9.output
    output:
        "data/processed-data/R9/events-train-norm.tsv"
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    log: 
        "results/logs/uncalled4-norm.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=500000,
        partition="main,zki,zki_mig", 
    shell:
        """
        python workflow/scripts/standardize-events.py {input} {output}
        """


rule eventalign_valid_norm_R9:
    input:
        rules.uncalled4_eventalign_valid_R9.output
    output:
        "data/processed-data/R9/events-valid-norm.tsv"
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    log: 
        "results/logs/uncalled4-norm-valid.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(7200),
        mem_mb=500000,
        partition="main,zki,zki_mig", 
    shell:
        """
        python workflow/scripts/standardize-events.py {input} {output}
        """



rule preprocess_train_seq2squiggle_R9:
    input:
        train=rules.eventalign_train_norm_R9.output,
        config="config/seq2squiggle-config-R9.yml",
    output:
        directory("results/train-human-R9/seq2squiggle-train"),
    conda:
        "../../envs/seq2squiggle-dev.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_preprocess_train_r9.txt"
    threads: 128
    log:
        "results/logs/preprocess_seq2squiggle.log",
    resources:
        disk_mb=100000,
        runtime=add_slack(3600),
        mem_mb=350000,
        partition="main"
    shell:
        """
        poetry -C resources/seq2squiggle/ install
        poetry -C resources/seq2squiggle/ run seq2squiggle preprocess {input.train} {output} --config {input.config}
        """
        # --max-chunks 100000000 are ~ 3mil chunks if mapping is 32 to 500


rule preprocess_valid_seq2squiggle_R9:
    input:
        valid=rules.eventalign_valid_norm_R9.output,
        config="config/seq2squiggle-config-R9.yml",
    output:
        directory("results/train-human-R9/seq2squiggle-valid"),
    conda:
        "../../envs/seq2squiggle-dev.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_preprocess_valid_r9.txt"
    threads: 128
    log:
        "results/logs/preprocess_seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(3600),
        mem_mb=100000,
        partition="main"
    shell:
        """
        poetry -C resources/seq2squiggle/ install
        poetry -C resources/seq2squiggle/ run seq2squiggle preprocess {input.valid} {output} --config {input.config}
        """


rule train_seq2squiggle_R9:
    input:
        train_dir=rules.preprocess_train_seq2squiggle_R9.output,
        valid_dir=rules.preprocess_valid_seq2squiggle_R9.output,
        config="config/seq2squiggle-config-R9.yml",
    output:
        "results/train-human-R9/seq2squiggle-logs/R9.ckpt",
    conda:
        "../../envs/seq2squiggle-dev.yml"
    benchmark:
        "results/benchmarks/seq2squiggle_training_r9.txt"
    threads: 64
    log:
        "results/logs/train_seq2squiggle.log",
    resources:
        disk_mb=50000, 
        runtime=add_slack(9000),
        mem_mb=150000, #300000
        slurm="gpus=8, nodes=1",
        partition="zki",     
    shell:
        """
        export HTTPS_PROXY="http://fw-bln.rki.local:8020"
        export https_proxy="http://fw-bln.rki.local:8020"
        poetry -C resources/seq2squiggle/ install
        poetry -C resources/seq2squiggle/ run seq2squiggle train {input.train_dir} {input.valid_dir} \
         --config {input.config} --model {output}
        """


rule test_bam_to_fastq_R9:
    input:
        rules.samtools_split_chr_uncalled4_R9.output.test,
    output:
        "data/processed-data/R9/reads_test.fastq"
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/bam_to_fastq.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        samtools bam2fq {input} > {output} 
        """


rule subsample_slow5_data_R9:
    input:
        fastq=rules.test_bam_to_fastq_R9.output,
        slow5=config["human_signal_R9"],
    output:
        read_l="data/processed-data/R9/reads_subset.txt",
        slow5="data/processed-data/R9/reads_subset.slow5",
    params:
        sample_rate=config["subsample-human-R9"],
    threads: 128
    log:
        "results/logs/subsample-slow5-human.log",
    resources:
        disk_mb=50000,
        mem_mb=100000,
    shell:
        # get read_list with x random headers from fastq and use slow5tools to extract these
        # removes @ symbol since this causes problems
        # awk '/^@/{{sub(/^@/, ""); print $1}}'
        """
        awk 'NR % 4 == 1 {{sub(/^@/, ""); print $1}}' {input.fastq} | shuf -n {params.sample_rate} > {output.read_l}
        resources/slow5tools-v1.1.0/slow5tools get -t {threads} {input.slow5} -l {output.read_l} -o {output.slow5} --skip
        resources/slow5tools-v1.1.0/slow5tools index {output.slow5}
        """

rule blow5_to_pod5_human_subset_R9:
    input:
        rules.subsample_slow5_data_R9.output.slow5,
    output:
        "data/processed-data/R9/reads_subset.pod5"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """

rule basecall_human_subset_R9:
    input:
        rules.blow5_to_pod5_human_subset_R9.output,
    output:
        "data/processed-data/R9/reads_subsample_pass.fastq",
    params:
        model=config["basecalling_model_R9"]
    threads: 128
    log:
        "results/logs/basecall_human_subset.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,
        slurm="gpus=2",
        partition="zki,main",   
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """


rule remove_header_fastq_subset_R9:
    input:
        rules.basecall_human_subset_R9.output,
    output:
        "data/processed-data/R9/reads_subsample_pass_fixed.fastq",
    threads: 64
    log:
        "results/logs/fix-fastq-header.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule remove_short_reads_subset_R9:
    input:
        rules.remove_header_fastq_subset_R9.output,
    output:
        "data/processed-data/R9/reads_subsample_pass_fixed_2.fastq",
    params:
        max_len=config["max_seq_len"]
    threads: 64
    log:
        "results/logs/fix-fastq-header.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{FS = "\\t" ; OFS = "\\n"}} {{header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= {params.max_len}) {{print header, seq, qheader, qseq}}}}' < {input} > {output}
        """

rule align_reference_subset_R9:
    input:
        seqs=rules.remove_short_reads_subset_R9.output,
        ref=config["human_ref"],
    output:
        "data/processed-data/R9/reads_subsample_pass.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align-reference.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """
