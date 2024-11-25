rule sam_stats_squigulator_ideal_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_squigulator_ideal_fly.output
    output:
        "results/d-melanogaster-noise/squigulator_ideal.npz",
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


rule sam_stats_squigulator_idealtime_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_squigulator_idealtime_fly.output
    output:
        "results/d-melanogaster-noise/squigulator_idealtime.npz",
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

rule sam_stats_squigulator_idealamp_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_squigulator_idealamp_fly.output
    output:
        "results/d-melanogaster-noise/squigulator_idealamp.npz",
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

rule sam_stats_squigulator_default_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_squigulator_default_fly.output
    output:
        "results/d-melanogaster-noise/squigulator_default.npz",
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


rule sam_stats_seq2squiggle_ideal_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_nonoise_fly.output
    output:
        "results/d-melanogaster-noise/seq2squiggle-no-noise.npz",
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

rule sam_stats_seq2squiggle_manualnoise_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-static-NS.npz",
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

rule sam_stats_seq2squiggle_noisesampler_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-NS.npz",
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


rule sam_stats_seq2squiggle_lengthregulator_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_lengthregulator_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-DS.npz",
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

rule sam_stats_seq2squiggle_lengthregulator_and_manualnoise_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_lengthregulator_and_manualnoise_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-DS-static-NS.npz",
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

rule sam_stats_seq2squiggle_manuallengthsampling_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_manuallengthsampling_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-static-DS.npz",
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

rule sam_stats_seq2squiggle_lengthregulator_and_noisesampelr_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_lengthregulator_and_noisesampler_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-DS-NS.npz",
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

rule sam_stats_seq2squiggle_mLR_mNS_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_seq2squiggle_mLR_mNS_fly.output,
    output:
        "results/d-melanogaster-noise/seq2squiggle-static-DS-static-NS.npz",
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


rule sam_stats_experimental_noise_fly:
    input:
        ref=config["d_melanogaster_ref"],
        sam=rules.align_reference_subset_fly.output
    output:
        "results/d-melanogaster-noise/experimental.npz",
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

rule plot_sam_stats_noise_fly_squigulator:
    input:
        squigulator_ideal=rules.sam_stats_squigulator_ideal_fly.output,
        squigulator_idealtime=rules.sam_stats_squigulator_idealtime_fly.output,
        squigulator_idealamp=rules.sam_stats_squigulator_idealamp_fly.output,
        squigulator_default=rules.sam_stats_squigulator_default_fly.output,
        experimental=rules.sam_stats_experimental_noise_fly.output,
    output:
        directory("results/d-melanogaster-noise/plots-squigulator"),
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
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.experimental} {input.squigulator_ideal} {input.squigulator_idealtime} {input.squigulator_idealamp} {input.squigulator_default} 
        """


rule plot_sam_stats_noise_fly_seq2squiggle:
    input:
        seq2squiggle_nonoise=rules.sam_stats_seq2squiggle_ideal_fly.output,
        seq2squiggle_manualnoise=rules.sam_stats_seq2squiggle_manualnoise_fly.output,
        seq2squiggle_noisesampler=rules.sam_stats_seq2squiggle_noisesampler_fly.output,
        seq2squiggle_lengthregulator=rules.sam_stats_seq2squiggle_lengthregulator_fly.output,
        seq2squiggle_lengthregulator_manualnoise=rules.sam_stats_seq2squiggle_lengthregulator_and_manualnoise_fly.output,
        seq2squiggle_manuallengthsampling=rules.sam_stats_seq2squiggle_manuallengthsampling_fly.output,
        seq2squiggle_LR_NS=rules.sam_stats_seq2squiggle_lengthregulator_and_noisesampelr_fly.output,
        seq2squiggle_mLR_mNS=rules.sam_stats_seq2squiggle_mLR_mNS_fly.output,
        experimental=rules.sam_stats_experimental_noise_fly.output,
    output:
        directory("results/d-melanogaster-noise/plots-seq2squiggle"),
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
        python workflow/scripts/basecalling-stats.py plot --outdir {output} {input.experimental} {input.seq2squiggle_nonoise} {input.seq2squiggle_manualnoise} {input.seq2squiggle_noisesampler} {input.seq2squiggle_lengthregulator} {input.seq2squiggle_lengthregulator_manualnoise} {input.seq2squiggle_manuallengthsampling} {input.seq2squiggle_LR_NS} {input.seq2squiggle_mLR_mNS}
        """

rule publish_results_noise_fly:
    input:
        plots_squigulator=rules.plot_sam_stats_noise_fly_squigulator.output,
        plots_seq2squiggle=rules.plot_sam_stats_noise_fly_seq2squiggle.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config.yml"
    output:
        directory("summary/d-melanogaster-noise/")
    log:
        "results/logs/publish-summary-noise.log",
    shell:
        # Save the config file as well since we often run the different parts
        # of the benchmark at different times
        """
        mkdir -p {output}
        cp -r {input.plots_squigulator} {output}
        cp -r {input.plots_seq2squiggle} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """