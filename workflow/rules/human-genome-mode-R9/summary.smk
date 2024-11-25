rule runtime_plot_gm_R9:
    input:
        squigulator=rules.run_squigulator_gm_R9.benchmark,
        seq2squiggle=rules.run_seq2squiggle_gm_R9.benchmark,
    output:
        "results/runtime-plot-R9.png"
    conda:
        "../../envs/plot-runtime.yml"
    shell:
        """
        python workflow/scripts/runtime_plot.py {input.squigulator} {input.seq2squiggle} {output}
        """


rule sam_stats_squigulator_gm_R9:
    input:
        ref=config["human_ref"],
        sam=rules.align_squigulator_gm_R9.output
    output:
        "results/human_gm_R9/squigulator.npz",
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

rule sam_stats_seq2squiggle_gm_R9:
    input:
        ref=config["human_ref"],
        sam=rules.align_seq2squiggle_gm_R9.output
    output:
        "results/human_gm_R9/seq2squiggle.npz",
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

rule sam_stats_experimental_gm_R9:
    input:
        ref=config["human_ref"],
        sam=rules.align_reference_subset_R9.output
    output:
        "results/human_gm_R9/experimental.npz",
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

rule sam_stats_deepsimulator_CI_gm_R9:
    input:
        ref=config["human_ref"],
        sam=rules.align_deepsimulator_CI_gm_R9.output
    output:
        "results/human_gm_R9/deepsimulator_CI.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-deepsimulator_CI.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule sam_stats_deepsimulator_CD_gm_R9:
    input:
        ref=config["human_ref"],
        sam=rules.align_deepsimulator_CD_gm_R9.output
    output:
        "results/human_gm_R9/deepsimulator_CD.npz",
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-deepsimulator_CD.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate --sam {input.sam} --ref {input.ref} --out {output}
        """

rule plot_sam_stats_gm_R9:
    input:
        seq2squiggle=rules.sam_stats_seq2squiggle_gm_R9.output,
        squigulator=rules.sam_stats_squigulator_gm_R9.output,
        experimental=rules.sam_stats_experimental_gm_R9.output,
        deepsimulator_CI=rules.sam_stats_deepsimulator_CI_gm_R9.output,
        deepsimulator_CD=rules.sam_stats_deepsimulator_CD_gm_R9.output,
    output:
        directory("results/human_gm_R9/plots"),
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
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.seq2squiggle} {input.squigulator} {input.experimental} {input.deepsimulator_CI} {input.deepsimulator_CD}
        """

rule publish_results_gm_R9:
    input:
        plots=rules.plot_sam_stats_gm_R9.output,
        runtime_plot=rules.runtime_plot_gm_R9.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config-R9.yml"
    output:
        directory("summary/human-genome-mode-R9/")
    log:
        "results/logs/publish-summary-gm.log",
    shell:
        # Save the config file as well since we often run the different parts
        # of the benchmark at different times
        """
        mkdir -p {output}
        cp -r {input.runtime_plot} {output}
        cp -r {input.plots} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """