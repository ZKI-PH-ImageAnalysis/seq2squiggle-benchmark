rule evaluate_before_finetuning:
    input:
        sam=rules.align_before_finetuning.output,
        ref=config["human_ref"],
    output:
        "results/fine-tuning-bonito/before-finetuning.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule evaluate_after_finetuning:
    input:
        sam=rules.align_after_finetuning.output
        ref=config["human_ref"],
    output:
        "results/fine-tuning-bonito/after-finetuning.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule plot_finetuning:
    input:
        before=rules.sam_stats_seq2squiggle_gm.output,
        after=rules.sam_stats_squigulator_gm.output,
    output:
        directory("results/human_gm/plots"),
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-plot.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.before} {input.after}
        """
