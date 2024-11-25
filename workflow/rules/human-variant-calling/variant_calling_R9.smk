rule run_squigulator_variant_R9:
    input:
        rules.generate_variant_genome.output.out,
    output:
        "results/variants_R9/squigulator_reads.slow5",
    params:
        sample_rate=config["subsample-variants-R9"],
        commands=config["variant_calling_R9"]["squigulator"],
    threads: 128
    log:
        "results/logs/run_squigulator_variant.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=150000,   
    shell:
        """
        resources/squigulator-v0.4.0/squigulator {params.commands} -n {params.sample_rate} {input} -o {output} 
        """


rule blow5_to_pod5_squigulator_variant_R9:
    input:
        rules.run_squigulator_variant_R9.output,
    output:
        "results/variants_R9/squigulator_reads.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """        

rule basecall_squigulator_variant_R9:
    input:
        rules.blow5_to_pod5_squigulator_variant_R9.output,
    output:
        "results/variants_R9/squigulator_reads.fastq"
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
        partition="zki",   
    shell:
        """
        ./resources/dorado-0.8.0-linux-x64/bin/dorado basecaller --emit-fastq {params.model} {input} > {output}
        """

rule remove_header_squigulator_fastq_variant_R9:
    input:
        rules.basecall_squigulator_variant_R9.output,
    output:
        "results/variants_R9/squigulator_reads_fixed.fastq",
    threads: 64
    log:
        "results/logs/remove_header_squigulator_fastq_variant.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        awk 'BEGIN {{ FS = " "; ORS = "\\n" }} /^@/ {{ print $1 }} !/^@/ {{ print }}' {input} > {output}
        """

rule align_squigulator_variant_R9:
    input:
        seqs=rules.remove_header_squigulator_fastq_variant_R9.output,
        ref=rules.extract_chr_from_genome.output,
    output:
        "results/variants_R9/squigulator_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align_squigulator_variant.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """


rule sort_index_squigulator_R9:
    input:
        rules.align_squigulator_variant_R9.output,
    output:
        "results/variants_R9/squigulator_reads.bam",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule variantcall_squigulator_R9:
    input:
        ref_genome=rules.extract_chr_from_genome.output,
        bam=rules.sort_index_squigulator_R9.output,
        model=config["clair3_model_path_R9"]
    output:
        directory("results/variants_R9/squigulator-vcf/"),
    params:
        commands=config["clair3_commands"]
    conda:
        "../../envs/clair3.yml"
    threads: 128
    log:
        "results/logs/variantcall_squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=300000,    
    shell:
        """
        run_clair3.sh --threads={threads} --bam_fn={input.bam} --ref_fn={input.ref_genome} --platform=ont \
        --model_path={input.model} --output={output} --sample_name=reads \
        {params.commands}
        """


rule run_seq2squiggle_variants_R9:
    input:
        model=rules.train_seq2squiggle_R9.output,
        fastq=rules.generate_variant_genome.output.out,
        config="config/seq2squiggle-config-R9.yml",
        installed=rules.install_poetry.output,
    output:
        "results/variants_R9/seq2squiggle_reads.slow5"
    params:
        commands=config["variant_calling_R9"]["seq2squiggle"],
        sample_rate=config["subsample-variants-R9"],
    conda:
        "../../envs/seq2squiggle-dev-copy.yml"
    log:
        "results/logs/run_seq2squiggle_variants_r9.log"
    threads: 128
    resources:
        disk_mb=500000,
        runtime=add_slack(1000),
        mem_mb=500000,
        slurm="gpus=1",
        partition="zki", 
    shell:
        """
        poetry -C resources/seq2squiggle/ run seq2squiggle predict --config {input.config} --model {input.model} {params.commands} -n {params.sample_rate} -o {output} {input.fastq}
        """

rule blow5_to_pod5_seq2squiggle_variant_R9:
    input:
        rules.run_seq2squiggle_variants_R9.output,
    output:
        "results/variants_R9/seq2squiggle_reads.pod5"
    conda:
        "../../envs/blue-crab.yml"
    threads: 128
    conda:
        "../../envs/minimap2.yml"
    shell:
        """
        blue-crab s2p {input} -o {output}
        """        

rule basecall_seq2squiggle_variants_R9:
    input:
        rules.blow5_to_pod5_seq2squiggle_variant_R9.output
    output:
        "results/variants_R9/seq2squiggle_reads.fastq"
    params:
        model=config["basecalling_model_R9"]
    threads: 32
    log:
        "results/logs/basecall_seq2squiggle_variants.log"
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

rule remove_header_seq2squiggle_fastq_variants_R9:
    input:
        rules.basecall_seq2squiggle_variants_R9.output,
    output:
        "results/variants_R9/seq2squiggle_reads_fixed.fastq",
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


rule align_seq2squiggle_variants_R9:
    input:
        seqs=rules.remove_header_seq2squiggle_fastq_variants_R9.output,
        ref=rules.extract_chr_from_genome.output,
    output:
        "results/variants_R9/seq2squiggle_reads.sam",
    params:
        commands=config["minimap2_commands"]
    conda:
        "../../envs/minimap2.yml"
    threads: 64
    log:
        "results/logs/align_seq2squiggle_variants.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        minimap2 {params.commands} -t {threads} {input.ref} {input.seqs} > {output}
        """


rule sort_index_seq2squiggle_R9:
    input:
        rules.align_seq2squiggle_variants_R9.output,
    output:
        "results/variants_R9/seq2squiggle_reads.bam",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule variantcall_seq2squiggle_R9:
    input:
        ref_genome=rules.extract_chr_from_genome.output,
        bam=rules.sort_index_seq2squiggle_R9.output,
        model=config["clair3_model_path_R9"]
    output:
        directory("results/variants_R9/seq2squiggle-vcf/"),
    params:
        commands=config["clair3_commands"]
    conda:
        "../../envs/clair3.yml"
    threads: 128
    log:
        "results/logs/variantcall_seq2squiggle_R9.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=200000,    
    shell:
        """
        run_clair3.sh --threads={threads} --bam_fn={input.bam} --ref_fn={input.ref_genome} --platform=ont \
        --model_path={input.model} --output={output} --sample_name=reads \
        {params.commands} 
        """



# DS CI
rule run_deepsimulator_CI_variants_R9:
    input:
        rules.generate_variant_genome.output.out,
    output:
        directory("results/variants_R9/deepsimulator_CI_reads/raw_reads"),
    conda:
        "../../envs/deepsimulator.yml"
    params:
        sample_rate=int(int(config["subsample-variants-R9"])/2), # necessary since DS is not able to process multi-fasta files
        run_params=config["variant_calling_R9"]["deepsimulator_CI"]
    threads: 64
    benchmark:
        "results/benchmarks/deepsimulator_CI_variants_R9.txt"
    log:
        "results/logs/deepsimulator_CI_variants-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=100000,   
    shell:
        """
        ./resources/DeepSimulator_benchmark/deep_simulator.sh -i {input} -o {output} -n {params.sample_rate} {params.run_params} -c {threads}
        """


rule fast5_to_blow5_CI_variants_R9:
    input:
        rules.run_deepsimulator_CI_variants_R9.output,
    output:
        directory("results/variants_R9/deepsimulator_CI_reads/blow5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        ./resources/slow5tools-v1.1.0/slow5tools f2s {input}/fast5/ -d {output}
        """         

rule blow5_to_pod5_CI_variants_R9:
    input:
        rules.fast5_to_blow5_CI_variants_R9.output,
    output:
        directory("results/variants_R9/deepsimulator_CI_reads/pod5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -d {output}
        """  


rule basecall_deepsimulator_CI_variants_R9:
    input:
        rules.blow5_to_pod5_CI_variants_R9.output,
    output:
        "results/variants_R9/deepsimulator_CI_reads.fastq"
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



rule align_deepsimulator_CI_variants_R9:
    input:
        seqs=rules.basecall_deepsimulator_CI_variants_R9.output,
        ref=rules.extract_chr_from_genome.output,
    output:
        "results/variants_R9/deepsimulator_CI_reads.sam",
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

rule sort_index_deepsimulator_CI_variants_R9:
    input:
        rules.align_deepsimulator_CI_variants_R9.output,
    output:
        "results/variants_R9/deepsimulator_CI_reads.bam",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_DS.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule variantcall_deepsimulator_CI_R9:
    input:
        ref_genome=rules.extract_chr_from_genome.output,
        bam=rules.sort_index_deepsimulator_CI_variants_R9.output,
        model=config["clair3_model_path_R9"]
    output:
        directory("results/variants_R9/deepsimulator-ci-vcf/"),
    params:
        commands=config["clair3_commands"]
    conda:
        "../../envs/clair3.yml"
    threads: 128
    log:
        "results/logs/variantcall_deepsimulatorci_R9.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=200000,    
    shell:
        """
        run_clair3.sh --threads={threads} --bam_fn={input.bam} --ref_fn={input.ref_genome} --platform=ont \
        --model_path={input.model} --output={output} --sample_name=reads \
        {params.commands} 
        """


# DS CD
rule run_deepsimulator_CD_variants_R9:
    input:
        rules.generate_variant_genome.output.out,
    output:
        directory("results/variants_R9/deepsimulator_CD_reads/raw_reads"),
    conda:
        "../../envs/deepsimulator.yml"
    params:
        sample_rate=int(int(config["subsample-variants-R9"])/2), # necessary since DS is not able to process multi-fasta files
        run_params=config["variant_calling_R9"]["deepsimulator_CD"]
    threads: 128
    benchmark:
        "results/benchmarks/deepsimulator_CD_variants_R9.txt"
    log:
        "results/logs/deepsimulator_CD_variants-simulation.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,   
    shell:
        """
        ./resources/DeepSimulator_benchmark/deep_simulator.sh -i {input} -o {output} -n {params.sample_rate} {params.run_params} -c {threads}
        """


rule fast5_to_blow5_CD_variants_R9:
    input:
        rules.run_deepsimulator_CD_variants_R9.output,
    output:
        directory("results/variants_R9/deepsimulator_CD_reads/blow5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        ./resources/slow5tools-v1.1.0/slow5tools f2s {input}/fast5/ -d {output}
        """         

rule blow5_to_pod5_CD_variants_R9:
    input:
        rules.fast5_to_blow5_CD_variants_R9.output,
    output:
        directory("results/variants_R9/deepsimulator_CD_reads/pod5")
    conda:
        "../../envs/seq2squiggle-dev.yml"
    threads: 128
    shell:
        """
        blue-crab s2p {input} -d {output}
        """  


rule basecall_deepsimulator_CD_variants_R9:
    input:
        rules.blow5_to_pod5_CD_variants_R9.output,
    output:
        "results/variants_R9/deepsimulator_CD_reads.fastq"
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



rule align_deepsimulator_CD_variants_R9:
    input:
        seqs=rules.basecall_deepsimulator_CD_variants_R9.output,
        ref=rules.extract_chr_from_genome.output,
    output:
        "results/variants_R9/deepsimulator_CD_reads.sam",
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

rule sort_index_deepsimulator_CD_variants_R9:
    input:
        rules.align_deepsimulator_CD_variants_R9.output,
    output:
        "results/variants_R9/deepsimulator_CD_reads.bam",
    conda:
        "../../envs/minimap2.yml"
    threads: 128
    log:
        "results/logs/sort_index_DS.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 
        samtools index {output}
        """

rule variantcall_deepsimulator_CD_R9:
    input:
        ref_genome=rules.extract_chr_from_genome.output,
        bam=rules.sort_index_deepsimulator_CD_variants_R9.output,
        model=config["clair3_model_path_R9"]
    output:
        directory("results/variants_R9/deepsimulator-cd-vcf/"),
    params:
        commands=config["clair3_commands"]
    conda:
        "../../envs/clair3.yml"
    threads: 128
    log:
        "results/logs/variantcall_deepsimulatorci_R9.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=200000,    
    shell:
        """
        run_clair3.sh --threads={threads} --bam_fn={input.bam} --ref_fn={input.ref_genome} --platform=ont \
        --model_path={input.model} --output={output} --sample_name=reads \
        {params.commands} 
        """



rule filter_vcf_files_R9:
    input:
        squigulator=rules.variantcall_squigulator_R9.output,
        seq2squiggle=rules.variantcall_seq2squiggle_R9.output,
        deepsimulator_ci=rules.variantcall_deepsimulator_CI_R9.output,
        deepsimulator_cd=rules.variantcall_deepsimulator_CD_R9.output,
    output:
        squigulator="results/variants_R9/squigulator-filtered.vcf.gz",
        seq2squiggle="results/variants_R9/seq2squiggle-filtered.vcf.gz",
        deepsimulator_ci="results/variants_R9/deepsimulator_ci-filtered.vcf.gz",
        deepsimulator_cd="results/variants_R9/deepsimulator_cd-filtered.vcf.gz",
    conda:
        "../../envs/clair3.yml"
    threads: 4
    log:
        "results/logs/rtg_eval_vcf.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        bcftools filter -i 'QUAL>0 && DP>1' {input.seq2squiggle}/merge_output.vcf.gz -Ov -o {output.seq2squiggle}
        tabix -f {output.seq2squiggle}
        bcftools filter -i 'QUAL>0 && DP>1' {input.squigulator}/merge_output.vcf.gz -Ov -o {output.squigulator}
        tabix -f {output.squigulator}
        bcftools filter -i 'QUAL>0 && DP>1' {input.deepsimulator_ci}/merge_output.vcf.gz -Ov -o {output.deepsimulator_ci}
        tabix -f {output.deepsimulator_ci}
        bcftools filter -i 'QUAL>0 && DP>1' {input.deepsimulator_cd}/merge_output.vcf.gz -Ov -o {output.deepsimulator_cd}
        tabix -f {output.deepsimulator_cd}
        """
    

rule rtg_eval_vcf_squigulator_R9:
    input:
        ref_genome=rules.rtg_index_format.output,
        squigulator=rules.filter_vcf_files_R9.output.squigulator,
        true=rules.filter_vcf_reference.output.vcf,
        bed=rules.filter_vcf_reference.output.bed,
    output:
        directory("results/variants_R9/eval-squigulator"),
    params:
        region=config["region_variants"]
    conda:
        "../../envs/clair3.yml"
    threads: 4
    log:
        "results/logs/rtg_eval_vcf.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        rtg vcfeval -b {input.true} -c {input.squigulator} \
         -t {input.ref_genome} --region {params.region} -e {input.bed} -f QUAL -o {output}
        """

rule rtg_eval_vcf_seq2squiggle_R9:
    input:
        ref_genome=rules.rtg_index_format.output,
        seq2squiggle=rules.filter_vcf_files_R9.output.seq2squiggle,
        true=rules.filter_vcf_reference.output.vcf,
        bed=rules.filter_vcf_reference.output.bed,
    output:
        directory("results/variants_R9/eval-seq2squiggle"),
    params:
        region=config["region_variants"]
    conda:
        "../../envs/clair3.yml"
    threads: 4
    log:
        "results/logs/rtg_eval_vcf_seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        rtg vcfeval -b {input.true} -c {input.seq2squiggle} \
         -t {input.ref_genome} --region {params.region} -e {input.bed} -f QUAL -o {output}
        """


rule rtg_eval_vcf_deepsimulator_ci_R9:
    input:
        ref_genome=rules.rtg_index_format.output,
        deepsimulator_ci=rules.filter_vcf_files_R9.output.deepsimulator_ci,
        true=rules.filter_vcf_reference.output.vcf,
        bed=rules.filter_vcf_reference.output.bed,
    output:
        directory("results/variants_R9/eval-deepsimulator-ci"),
    params:
        region=config["region_variants"]
    conda:
        "../../envs/clair3.yml"
    threads: 4
    log:
        "results/logs/rtg_eval_vcf.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        rtg vcfeval -b {input.true} -c {input.deepsimulator_ci} \
         -t {input.ref_genome} --region {params.region} -e {input.bed} -f QUAL -o {output}
        """

rule rtg_eval_vcf_deepsimulator_cd_R9:
    input:
        ref_genome=rules.rtg_index_format.output,
        deepsimulator_cd=rules.filter_vcf_files_R9.output.deepsimulator_cd,
        true=rules.filter_vcf_reference.output.vcf,
        bed=rules.filter_vcf_reference.output.bed,
    output:
        directory("results/variants_R9/eval-deepsimulator-cd"),
    params:
        region=config["region_variants"]
    conda:
        "../../envs/clair3.yml"
    threads: 4
    log:
        "results/logs/rtg_eval_vcf.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        rtg vcfeval -b {input.true} -c {input.deepsimulator_cd} \
         -t {input.ref_genome} --region {params.region} -e {input.bed} -f QUAL -o {output}
        """


rule annotate_FP_variants_squigulator_R9:
    input:
        squigulator=rules.rtg_eval_vcf_squigulator_R9.output,
        homopolymer_regions=rules.extract_homopolymer_regions.output,
        repeat_regions=rules.extract_repeat_regions.output.bed,
    output:
        "results/variants_R9/regions_squigulator.txt"
    conda:
        "../../envs/bcftools.yml"
    threads: 16
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        handle_headers() {{
            input_file=$1
            output_file=$2
            original_vcf=$3

            if ! grep -q '^#' "$input_file"; then
                bcftools view -h "$original_vcf" > results/variants/header.txt
                cat results/variants/header.txt "$input_file" > "$output_file"
            else
                mv "$input_file" "$output_file"
            fi
        }}

        # Process false positives
        bedtools intersect -a {input.squigulator}/fp.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_FP_squigulator_homop.vcf
        bedtools intersect -a {input.squigulator}/fp.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_FP_squigulator_repeat.vcf

        handle_headers results/variants/temp_FP_squigulator_homop.vcf results/variants/FP_squigulator_homop.vcf {input.squigulator}/fp.vcf.gz
        handle_headers results/variants/temp_FP_squigulator_repeat.vcf results/variants/FP_squigulator_repeat.vcf {input.squigulator}/fp.vcf.gz

        count_homopolymers=$(bcftools view -H results/variants/FP_squigulator_homop.vcf | wc -l)
        count_strs=$(bcftools view -H results/variants/FP_squigulator_repeat.vcf | wc -l)

        cat {input.homopolymer_regions} {input.repeat_regions} | sort -k1,1 -k2,2n | bedtools merge -i - > results/variants/combined_homopolymer_str_regions.bed

        bedtools intersect -a {input.squigulator}/fp.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_FP_squigulator_other.vcf
        handle_headers results/variants/temp_FP_squigulator_other.vcf results/variants/FP_squigulator_other.vcf {input.squigulator}/fp.vcf.gz
        count_other=$(bcftools view -H results/variants/FP_squigulator_other.vcf | wc -l)

        # Process false negatives
        bedtools intersect -a {input.squigulator}/fn.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_FN_squigulator_homop.vcf
        bedtools intersect -a {input.squigulator}/fn.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_FN_squigulator_repeat.vcf

        handle_headers results/variants/temp_FN_squigulator_homop.vcf results/variants/FN_squigulator_homop.vcf {input.squigulator}/fn.vcf.gz
        handle_headers results/variants/temp_FN_squigulator_repeat.vcf results/variants/FN_squigulator_repeat.vcf {input.squigulator}/fn.vcf.gz

        bedtools intersect -a {input.squigulator}/fn.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_FN_squigulator_other.vcf
        handle_headers results/variants/temp_FN_squigulator_other.vcf results/variants/FN_squigulator_other.vcf {input.squigulator}/fn.vcf.gz
        count_other_FN=$(bcftools view -H results/variants/FN_squigulator_other.vcf | wc -l)

        count_homopolymers_FN=$(bcftools view -H results/variants/FN_squigulator_homop.vcf | wc -l)
        count_strs_FN=$(bcftools view -H results/variants/FN_squigulator_repeat.vcf | wc -l)

        # Process True postives
        bedtools intersect -a {input.squigulator}/tp.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_TP_squigulator_homop.vcf
        bedtools intersect -a {input.squigulator}/tp.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_TP_squigulator_repeat.vcf

        handle_headers results/variants/temp_TP_squigulator_homop.vcf results/variants/TP_squigulator_homop.vcf {input.squigulator}/tp.vcf.gz
        handle_headers results/variants/temp_TP_squigulator_repeat.vcf results/variants/TP_squigulator_repeat.vcf {input.squigulator}/tp.vcf.gz

        bedtools intersect -a {input.squigulator}/tp.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_TP_squigulator_other.vcf
        handle_headers results/variants/temp_TP_squigulator_other.vcf results/variants/TP_squigulator_other.vcf {input.squigulator}/tp.vcf.gz
        
        count_other_TP=$(bcftools view -H results/variants/TP_squigulator_other.vcf | wc -l)
        count_homopolymers_TP=$(bcftools view -H results/variants/TP_squigulator_homop.vcf | wc -l)
        count_strs_TP=$(bcftools view -H results/variants/TP_squigulator_repeat.vcf | wc -l)


        # Output the results
        echo "False positives due to homopolymers: $count_homopolymers" > {output}
        echo "False positives due to STRs: $count_strs" >> {output}
        echo "False positives due to other regions: $count_other" >> {output}
        echo "False negatives due to homopolymers: $count_homopolymers_FN" >> {output}
        echo "False negatives due to STRs: $count_strs_FN" >> {output}
        echo "False negatives due to other regions: $count_other_FN" >> {output}

        echo "True positives due to homopolymers: $count_homopolymers_TP" >> {output}
        echo "True positives due to STRs: $count_strs_TP" >> {output}
        echo "True positives due to other regions: $count_other_TP" >> {output}

        rm -rf results/variants/temp*squigulator
        rm -rf results/variants/TP*squigulator
        rm -rf results/variants/FP*squigulator
        rm -rf results/variants/FN*squigulator
        """

rule annotate_FP_variants_seq2squiggle_R9:
    input:
        seq2squiggle=rules.rtg_eval_vcf_seq2squiggle_R9.output,
        homopolymer_regions=rules.extract_homopolymer_regions.output,
        repeat_regions=rules.extract_repeat_regions.output.bed,
    output:
        "results/variants_R9/regions_seq2squiggle.txt"
    conda:
        "../../envs/bcftools.yml"
    threads: 16
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000,    
    shell:
        """
        handle_headers() {{
            input_file=$1
            output_file=$2
            original_vcf=$3

            if ! grep -q '^#' "$input_file"; then
                bcftools view -h "$original_vcf" > results/variants/header.txt
                cat results/variants/header.txt "$input_file" > "$output_file"
            else
                mv "$input_file" "$output_file"
            fi
        }}

        # Process false positives
        bedtools intersect -a {input.seq2squiggle}/fp.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_FP_seq2squiggle_homop.vcf
        bedtools intersect -a {input.seq2squiggle}/fp.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_FP_seq2squiggle_repeat.vcf

        handle_headers results/variants/temp_FP_seq2squiggle_homop.vcf results/variants/FP_seq2squiggle_homop.vcf {input.seq2squiggle}/fp.vcf.gz
        handle_headers results/variants/temp_FP_seq2squiggle_repeat.vcf results/variants/FP_seq2squiggle_repeat.vcf {input.seq2squiggle}/fp.vcf.gz

        count_homopolymers=$(bcftools view -H results/variants/FP_seq2squiggle_homop.vcf | wc -l)
        count_strs=$(bcftools view -H results/variants/FP_seq2squiggle_repeat.vcf | wc -l)

        cat {input.homopolymer_regions} {input.repeat_regions} | sort -k1,1 -k2,2n | bedtools merge -i - > results/variants/combined_homopolymer_str_regions.bed

        bedtools intersect -a {input.seq2squiggle}/fp.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_FP_seq2squiggle_other.vcf
        handle_headers results/variants/temp_FP_seq2squiggle_other.vcf results/variants/FP_seq2squiggle_other.vcf {input.seq2squiggle}/fp.vcf.gz
        count_other=$(bcftools view -H results/variants/FP_seq2squiggle_other.vcf | wc -l)

        # Process false negatives
        bedtools intersect -a {input.seq2squiggle}/fn.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_FN_seq2squiggle_homop.vcf
        bedtools intersect -a {input.seq2squiggle}/fn.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_FN_seq2squiggle_repeat.vcf

        handle_headers results/variants/temp_FN_seq2squiggle_homop.vcf results/variants/FN_seq2squiggle_homop.vcf {input.seq2squiggle}/fn.vcf.gz
        handle_headers results/variants/temp_FN_seq2squiggle_repeat.vcf results/variants/FN_seq2squiggle_repeat.vcf {input.seq2squiggle}/fn.vcf.gz

        bedtools intersect -a {input.seq2squiggle}/fn.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_FN_seq2squiggle_other.vcf
        handle_headers results/variants/temp_FN_seq2squiggle_other.vcf results/variants/FN_seq2squiggle_other.vcf {input.seq2squiggle}/fn.vcf.gz
        count_other_FN=$(bcftools view -H results/variants/FN_seq2squiggle_other.vcf | wc -l)

        count_homopolymers_FN=$(bcftools view -H results/variants/FN_seq2squiggle_homop.vcf | wc -l)
        count_strs_FN=$(bcftools view -H results/variants/FN_seq2squiggle_repeat.vcf | wc -l)

        # Process True postives
        bedtools intersect -a {input.seq2squiggle}/tp.vcf.gz -b {input.homopolymer_regions} -wa -u > results/variants/temp_TP_seq2squiggle_homop.vcf
        bedtools intersect -a {input.seq2squiggle}/tp.vcf.gz -b {input.repeat_regions} -wa -u > results/variants/temp_TP_seq2squiggle_repeat.vcf

        handle_headers results/variants/temp_TP_seq2squiggle_homop.vcf results/variants/TP_seq2squiggle_homop.vcf {input.seq2squiggle}/tp.vcf.gz
        handle_headers results/variants/temp_TP_seq2squiggle_repeat.vcf results/variants/TP_seq2squiggle_repeat.vcf {input.seq2squiggle}/tp.vcf.gz

        bedtools intersect -a {input.seq2squiggle}/tp.vcf.gz -b results/variants/combined_homopolymer_str_regions.bed -v > results/variants/temp_TP_seq2squiggle_other.vcf
        handle_headers results/variants/temp_TP_seq2squiggle_other.vcf results/variants/TP_seq2squiggle_other.vcf {input.seq2squiggle}/tp.vcf.gz
        
        count_other_TP=$(bcftools view -H results/variants/TP_seq2squiggle_other.vcf | wc -l)
        count_homopolymers_TP=$(bcftools view -H results/variants/TP_seq2squiggle_homop.vcf | wc -l)
        count_strs_TP=$(bcftools view -H results/variants/TP_seq2squiggle_repeat.vcf | wc -l)


        # Output the results
        echo "False positives due to homopolymers: $count_homopolymers" > {output}
        echo "False positives due to STRs: $count_strs" >> {output}
        echo "False positives due to other regions: $count_other" >> {output}
        echo "False negatives due to homopolymers: $count_homopolymers_FN" >> {output}
        echo "False negatives due to STRs: $count_strs_FN" >> {output}
        echo "False negatives due to other regions: $count_other_FN" >> {output}

        echo "True positives due to homopolymers: $count_homopolymers_TP" >> {output}
        echo "True positives due to STRs: $count_strs_TP" >> {output}
        echo "True positives due to other regions: $count_other_TP" >> {output}

        rm -rf results/variants/temp*seq2squiggle
        rm -rf results/variants/TP*seq2squiggle
        rm -rf results/variants/FP*seq2squiggle
        rm -rf results/variants/FN*seq2squiggle
        """




rule plot_variants_region_R9:
    input:
        squigulator=rules.annotate_FP_variants_squigulator_R9.output,
        seq2squiggle=rules.annotate_FP_variants_seq2squiggle_R9.output,
    output:
        "results/variants_R9/SupplementaryFig09.png"
    conda:
        "../../envs/sam-stats.yml"
    threads: 16
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=100000, 
    shell:
        """
        python workflow/scripts/variants-region-plot.py {input.seq2squiggle} {input.squigulator} {output}
        """


rule sam_stats_seq2squiggle_variants_R9:
    input:
        ref=rules.extract_chr_from_genome.output,
        sam=rules.align_seq2squiggle_variants_R9.output
    output:
        "results/variants_R9/seq2squiggle.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule sam_stats_squigulator_variants_R9:
    input:
        ref=rules.extract_chr_from_genome.output,
        sam=rules.align_squigulator_variant_R9.output
    output:
        "results/variants_R9/squigulator.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=20000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule plot_sam_stats_variants_R9:
    input:
        seq2squiggle=rules.sam_stats_seq2squiggle_variants_R9.output,
        squigulator=rules.sam_stats_squigulator_variants_R9.output,
    output:
        directory("results/variants_R9/plots"),
    threads: 1
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-plot.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=20000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.seq2squiggle} {input.squigulator}
        """


rule plot_rtg_weighted_roc_R9:
    input:
        squigulator=rules.rtg_eval_vcf_squigulator_R9.output,
        seq2squiggle=rules.rtg_eval_vcf_seq2squiggle_R9.output,
        deepsimulator_ci=rules.rtg_eval_vcf_deepsimulator_ci_R9.output,
        deepsimulator_cd=rules.rtg_eval_vcf_deepsimulator_cd_R9.output,
    output:
        directory("results/variants_R9/weighted-plots")
    conda:
        "../../envs/sam-stats.yml"
    threads: 1
    log:
        "results/logs/sort_index_bam_valid.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=10000,    
    shell:
        """
        python workflow/scripts/rocplot.py --outdir {output} {input.seq2squiggle}/weighted_roc.tsv.gz \
        {input.squigulator}/weighted_roc.tsv.gz {input.deepsimulator_cd}/weighted_roc.tsv.gz \
        {input.deepsimulator_ci}/weighted_roc.tsv.gz \
        --tools seq2squiggle squigulator "DeepSimulator (context-dependent)" "DeepSimulator (context-independent)"
        """


rule plot_rtg_snp_roc_R9:
    input:
        squigulator=rules.rtg_eval_vcf_squigulator_R9.output,
        seq2squiggle=rules.rtg_eval_vcf_seq2squiggle_R9.output,
        deepsimulator_ci=rules.rtg_eval_vcf_deepsimulator_ci_R9.output,
        deepsimulator_cd=rules.rtg_eval_vcf_deepsimulator_cd_R9.output,
    output:
        directory("results/variants_R9/snp-plots")
    conda:
        "../../envs/sam-stats.yml"
    threads: 1
    log:
        "results/logs/sort_index_bam_valid.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=10000,    
    shell:
        """
        python workflow/scripts/rocplot.py --outdir {output} {input.seq2squiggle}/snp_roc.tsv.gz \
        {input.squigulator}/snp_roc.tsv.gz {input.deepsimulator_cd}/snp_roc.tsv.gz \
        {input.deepsimulator_ci}/snp_roc.tsv.gz \
        --tools seq2squiggle squigulator "DeepSimulator (context-dependent)" "DeepSimulator (context-independent)"
        """


rule plot_rtg_indel_roc_R9:
    input:
        squigulator=rules.rtg_eval_vcf_squigulator_R9.output,
        seq2squiggle=rules.rtg_eval_vcf_seq2squiggle_R9.output,
        deepsimulator_ci=rules.rtg_eval_vcf_deepsimulator_ci_R9.output,
        deepsimulator_cd=rules.rtg_eval_vcf_deepsimulator_cd_R9.output,
    output:
        directory("results/variants_R9/indel-plots")
    conda:
        "../../envs/sam-stats.yml"
    threads: 1
    log:
        "results/logs/sort_index_bam_valid.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=10000,    
    shell:
        """
        python workflow/scripts/rocplot.py --outdir {output} {input.seq2squiggle}/non_snp_roc.tsv.gz \
        {input.squigulator}/non_snp_roc.tsv.gz {input.deepsimulator_cd}/non_snp_roc.tsv.gz \
        {input.deepsimulator_ci}/non_snp_roc.tsv.gz \
        --tools seq2squiggle squigulator "DeepSimulator (context-dependent)" "DeepSimulator (context-independent)"
        """

rule publish_results_variants_R9:
    input:
        weighted_roc=rules.plot_rtg_weighted_roc_R9.output,
        snp_roc=rules.plot_rtg_snp_roc_R9.output,
        indel_roc=rules.plot_rtg_indel_roc_R9.output,
        plots=rules.plot_sam_stats_variants_R9.output,
        #regions_variants=rules.plot_variants_region_R9.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config-R9.yml"
    output:
        directory("summary/variants_R9/")
    log:
        "results/logs/publish_results_variants.log",
    shell:
        """
        mkdir -p {output}
        cp -r {input.weighted_roc} {output}
        cp -r {input.snp_roc} {output}
        cp -r {input.indel_roc} {output}
        cp -r {input.plots} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """
        # cp -r {input.regions_variants} {output}


rule plot_group_SNP:
    input:
        squigulator_R9=rules.rtg_eval_vcf_squigulator_R9.output,
        seq2squiggle_R9=rules.rtg_eval_vcf_seq2squiggle_R9.output,
        deepsimulator_ci_R9=rules.rtg_eval_vcf_deepsimulator_ci_R9.output,
        deepsimulator_cd_R9=rules.rtg_eval_vcf_deepsimulator_cd_R9.output,
        seq2squiggle_R10=rules.rtg_eval_vcf_seq2squiggle.output,
        squigulator_R10=rules.rtg_eval_vcf_squigulator.output,
    output:
        directory("summary/variants-all/snp-plots")
    conda:
        "../../envs/sam-stats.yml"
    threads: 1
    log:
        "results/logs/plot-group-snp.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(1800),
        mem_mb=10000,    
    shell:
        """
        python workflow/scripts/rocplot_groupplot.py --outdir {output} --r9_files {input.seq2squiggle_R9}/snp_roc.tsv.gz \
        {input.squigulator_R9}/snp_roc.tsv.gz {input.deepsimulator_cd_R9}/snp_roc.tsv.gz \
        {input.deepsimulator_ci_R9}/snp_roc.tsv.gz \
        --r9_tools seq2squiggle squigulator "DeepSimulator (context-dependent)" "DeepSimulator (context-independent)" \
        --r10_files {input.seq2squiggle_R10}/snp_roc.tsv.gz {input.squigulator_R10}/snp_roc.tsv.gz \
        --r10_tools seq2squiggle squigulator
        """

