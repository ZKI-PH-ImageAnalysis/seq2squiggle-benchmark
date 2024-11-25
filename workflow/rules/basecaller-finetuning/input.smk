rule basecall_bonito_before_finetuning:
    input:
        rules.blow5_to_pod5_human_subset.output,
    output:
        "results/fine-tuning-bonito/reads_default_bonnito.fastq",
    params:
        model=config["bonito_model"]
    conda:
        "../../envs/bonito.yml"
    threads: 64
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=50000,
        slurm="gpus=1",
        partition="zki,main",   
    shell:
        """
        bonito basecaller --emit-fastq {params.model} {input} > {output}
        """

rule align_before_finetuning:
    input:
        seqs=rules.basecall_bonito_before_finetuning.output,
        ref=config["human_ref"],
    output:
        "results/fine-tuning-bonito/reads_default.sam",
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


rule move_seq2squiggle_reads:
    input:
        # reads from seq2squiggle
    output:
        directory("results/fine-tuning-reads/all-reads")
    shell:
        """
        mdir {output}
        mv {input} {output}
        """

# TODO Do this for multiple noise modes of seq2squiggle and for multiple  genomes..
# Use maybe wildcard
rule basecall_for_finetune_bonito:
    input:
        seqs=rules.move_seq2squiggle_reads.output,
        ref=config["human_ref"],
    output:
        "results/fine-tuning-reads/ctc-data/basecalls.sam",
    params:
        model=config["basecalling_model"]
    conda:
        "../../envs/bonito.yml"
    threads: 64
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=50000,
        slurm="gpus=1",
        partition="zki,main",   
    shell:
        """
        bonito basecaller {params.model} --save-ctc --reference {input.ref} {input.seqs} > {output}
        """

rule finetune_bonito:
    input:
        "results/fine-tuning-bonito/ctc-data"
    output: 
        "results/fine-tuning-bonito/fine-tuned-model/"
    params:
        epochs=config["bonito_epochs"]
        lr=config["bonito_lr"]
        model=config["bonito_model"]
    conda:
        "../../envs/bonito.yml"
    threads: 64
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=50000,
        slurm="gpus=1",
        partition="zki,main",  
    shell:
        """
        bonito train --epochs {params.epochs} --lr {params.lr} --pretrained {params.model} --directory {input} {output}
        """

rule basecall_bonito_after_finetuning:
    input:
        seqs=rules.blow5_to_pod5_human_subset.output,
        model=rules.finetune_bonito.output,
    output:
        "results/fine-tuning-bonito/reads_default_bonnito.fastq",
    threads: 64
    conda:
        "../../envs/bonito.yml"
    resources:
        disk_mb=50000,
        runtime=add_slack(2100),
        mem_mb=50000,
        slurm="gpus=1",
        partition="zki,main",   
    shell:
        """
        bonito basecaller --emit-fastq {input.model} {input} > {output}
        """

rule align_after_finetuning:
    input:
        seqs=rules.basecall_bonito_after_finetuning.output,
        ref=config["human_ref"],
    output:
        "results/fine-tuning-bonito/reads_finetuned.sam",
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
