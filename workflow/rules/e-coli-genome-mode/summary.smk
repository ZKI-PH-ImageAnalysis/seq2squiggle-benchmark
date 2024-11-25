rule runtime_plot_gm_ecoli:
    input:
        squigulator=rules.run_squigulator_gm.benchmark,
        seq2squiggle=rules.run_seq2squiggle_gm.benchmark,
    output:
        "results/runtime-plot-ecoli.png"
    conda:
        "../../envs/plot-runtime.yml"
    shell:
        """
        python workflow/scripts/runtime_plot.py {input.squigulator} {input.seq2squiggle} {output}
        """


rule sam_stats_squigulator_gm_ecoli:
    input:
        ref=config["e_coli_ref"],
        sam=rules.align_squigulator_gm_ecoli.output
    output:
        "results/e-coli/squigulator.npz",
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

rule sam_stats_seq2squiggle_gm_ecoli:
    input:
        ref=config["e_coli_ref"],
        sam=rules.align_seq2squiggle_gm_ecoli.output
    output:
        "results/e-coli/seq2squiggle.npz",
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

rule sam_stats_experimental_gm_ecoli:
    input:
        ref=config["e_coli_ref"],
        sam=rules.align_reference_subset_ecoli.output
    output:
        "results/e-coli/experimental.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-experimental.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,   
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule plot_sam_stats_gm_ecoli:
    input:
        seq2squiggle=rules.sam_stats_seq2squiggle_gm_ecoli.output,
        squigulator=rules.sam_stats_squigulator_gm_ecoli.output,
        experimental=rules.sam_stats_experimental_gm_ecoli.output,
    output:
        directory("results/e-coli/plots"),
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
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.seq2squiggle} {input.squigulator} {input.experimental}
        """

rule publish_results_gm_ecoli:
    input:
        plots=rules.plot_sam_stats_gm_ecoli.output,
        runtime_plot=rules.runtime_plot_gm_ecoli.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config.yml"
    output:
        directory("summary/e-coli/")
    log:
        "results/logs/publish-summary-gm.log",
    shell:
        """
        mkdir -p {output}
        cp -r {input.runtime_plot} {output}
        cp -r {input.plots} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """