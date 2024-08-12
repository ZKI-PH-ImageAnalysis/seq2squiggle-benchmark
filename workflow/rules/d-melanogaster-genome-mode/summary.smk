rule runtime_plot_gm_fly:
    input:
        squigulator=rules.run_squigulator_gm.benchmark,
        seq2squiggle=rules.run_seq2squiggle_gm.benchmark,
    output:
        "results/runtime-plot-fly.png"
    conda:
        "../../envs/plot-runtime.yml"
    shell:
        """
        python workflow/scripts/runtime_plot.py {input.squigulator} {input.seq2squiggle} {output}
        """


rule sam_stats_squigulator_gm_fly:
    input:
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        sam=rules.align_squigulator_gm_fly.output
    output:
        "results/d-melanogaster/squigulator.npz",
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

rule sam_stats_seq2squiggle_gm_fly:
    input:
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        sam=rules.align_seq2squiggle_gm_fly.output
    output:
        "results/d-melanogaster/seq2squiggle.npz",
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

rule sam_stats_experimental_gm_fly:
    input:
        ref="data/d-melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
        sam=rules.align_reference_subset_fly.output
    output:
        "results/d-melanogaster/experimental.npz",
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

rule plot_sam_stats_gm_fly:
    input:
        seq2squiggle=rules.sam_stats_seq2squiggle_gm_fly.output,
        squigulator=rules.sam_stats_squigulator_gm_fly.output,
        experimental=rules.sam_stats_experimental_gm_fly.output,
    output:
        directory("results/d-melanogaster/plots"),
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

rule publish_results_gm_fly:
    input:
        plots=rules.plot_sam_stats_gm_fly.output,
        runtime_plot=rules.runtime_plot_gm_fly.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config.yml"
    output:
        directory("summary/d-melanogaster/")
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