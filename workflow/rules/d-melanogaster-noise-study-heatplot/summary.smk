rule sam_stats:
    input:
        ref=config["d_melanogaster_ref"],
        sam=lambda wildcards: f"results/dwell-{wildcards.dwell_std}_noise-{wildcards.noise_std}/{wildcards.tool}_reads.sam"
    output:
        "results/dwell-{dwell_std}_noise-{noise_std}/{tool}.npz"
    params:
        noise_std=lambda wildcards: wildcards.noise_std,
        dwell_std=lambda wildcards: wildcards.dwell_std,
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/dwell-{dwell_std}_noise-{noise_std}-{tool}.log"
    resources:
        disk_mb=50000,
        runtime=add_slack(900),
        mem_mb=200000
    shell:
        """
        python workflow/scripts/basecalling-stats.py evaluate \
            --sam {input.sam} --ref {input.ref} --out {output}
        """


rule plot_sam_stats_heatplot:
    input:
        # Gather all npz files for each tool
        lambda wildcards: sorted(set(expand("results/dwell-{dwell_std}_noise-{noise_std}/{tool}.npz", 
                                           dwell_std=[c[0] for c in combinations], 
                                           noise_std=[c[1] for c in combinations], 
                                           tool=[wildcards.tool])))
    output:
        "summary/{tool}_heatmap.png"  # One heatmap image per tool
    params:
        dwell_std_range=','.join(map(str, sorted(set([c[0] for c in combinations])))),
        noise_std_range=','.join(map(str, sorted(set([c[1] for c in combinations])))),
        colormap=lambda wildcards: "Oranges" if wildcards.tool == "squigulator" else "Blues",
        toolname=lambda wildcards: wildcards.tool  
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/plot-{tool}-heatmap.log"
    shell:
        """
        python workflow/scripts/basecalling-stats.py plot-heatmap \
            --out {output} \
            --dwell_stds {params.dwell_std_range} \
            --noise_stds {params.noise_std_range} \
            --colormap {params.colormap} \
            --toolname {params.toolname} \
            {input} 
        """