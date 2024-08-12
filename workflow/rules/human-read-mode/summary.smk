rule index_subsample_slow5:
    input:
        rules.subsample_slow5_data.output.slow5
    output:
        "data/zymo-human/PGXX22563_pcr/PGXX22563_reads_subsample.slow5.idx"
    threads: 128
    conda: 
        "../../envs/signal-sim.yml"
    log:
        "results/logs/index_slow5.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(400),
        mem_mb=200000,  
    shell:
        """
        slow5tools index {input}
        """


rule calc_dtw_squigulator:
    input:
        sim=rules.fix_squigulator_slow5.output,
        ref=rules.subsample_slow5_data.output.slow5,
        index=rules.index_subsample_slow5.output,
    output:
        "results/human/squigulator_dtw.tsv",
    threads: 128
    conda:
        "../../envs/signal-sim.yml"
    params:
        dtw_seqs=config["dtw_seqs"],
    log:
        "results/logs/calc-dtw-squigulator.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,   
    shell:
        """
        python workflow/scripts/signal-similarity.py compare --reference {input.ref} {input.sim} --out {output} --head_n {params.dtw_seqs} --norm True --dtw_impl "linmdtw-fastdtw"
        """


rule calc_dtw_seq2squiggle:
    input:
        sim=rules.run_seq2squiggle.output,
        ref=rules.subsample_slow5_data.output.slow5,
        index=rules.index_subsample_slow5.output,
    output:
        "results/human/seq2squiggle_dtw.tsv",
    threads: 128
    conda:
        "../../envs/signal-sim.yml"
    params:
        dtw_seqs=config["dtw_seqs"],
    log:
        "results/logs/calc-dtw-seq2squiggle.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,   
    shell:
        """
        python workflow/scripts/signal-similarity.py compare --reference {input.ref} {input.sim} --out {output} --head_n {params.dtw_seqs} --norm True --dtw_impl "linmdtw-fastdtw"
        """

rule plot_dtw:
    input:
        seq2squiggle=rules.calc_dtw_seq2squiggle.output,
        squigulator=rules.calc_dtw_squigulator.output
    output:
        "results/human/dtw-analysis.png",
    threads: 128
    conda:
        "../../envs/signal-sim.yml"
    log:
        "results/logs/plot-dtw.log",
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000,  
    shell:
        """
        python workflow/scripts/signal-similarity.py plot --seq2squiggle {input.seq2squiggle} --squigulator {input.squigulator} --out {output}
        """

rule sam_stats_squigulator:
    input:
        ref="data/zymo-human/GRCh38_no_alt.fna",
        sam=rules.align_squigulator.output
    output:
        "results/human/squigulator.npz",
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

rule sam_stats_seq2squiggle:
    input:
        ref="data/zymo-human/GRCh38_no_alt.fna",
        sam=rules.align_seq2squiggle.output
    output:
        "results/human/seq2squiggle.npz",
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

rule sam_stats_experimental:
    input:
        ref="data/zymo-human/GRCh38_no_alt.fna",
        sam=rules.align_reference_subset.output
    output:
        "results/human/experimental.npz",
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

rule plot_sam_stats:
    input:
        seq2squiggle=rules.sam_stats_seq2squiggle.output,
        squigulator=rules.sam_stats_squigulator.output,
        experimental=rules.sam_stats_experimental.output
    output:
        directory("results/human/plots"),
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

rule publish_results:
    input:
        plots=rules.plot_sam_stats.output,
        #f5c_plots=rules.analyse_f5c_file.output,
        dtw_plots=rules.plot_dtw.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config.yml"
    output:
        directory("summary/human-read-mode/")
    log:
        "results/logs/publish-summary.log",
    shell:
        """
        mkdir -p {output}
        cp -r {input.plots} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """
