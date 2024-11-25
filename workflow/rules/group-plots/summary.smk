rule run_seq2squiggle_Figure_01:
    input:
        model=rules.train_seq2squiggle.output,
        fasta="data/zymo-human/Figure01_seq2squiggle_input.fasta",
        config="config/seq2squiggle-config.yml",
    output:
        "results/Figure01/Figure01_prediction.slow5"
    conda:
        "../../envs/seq2squiggle.yml"
    log:
        "results/logs/run_seq2squiggle_Figure_01.log"
    threads: 8
    resources:
        disk_mb=50000,
        runtime=add_slack(100),
        mem_mb=35000,
        slurm="gpus=1",
        partition="zki", 
    shell:
        """
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --read-input --config {input.config} --model {input.model}  --ideal-event-length -1.0 --duration-sampling True --noise-sampling True --noise-std 1.5 -o {output} {input.fasta}
        """


rule run_seq2squiggle_Figure_01_ideal:
    input:
        model=rules.train_seq2squiggle.output,
        fasta="data/zymo-human/Figure01_seq2squiggle_input.fasta",
        config="config/seq2squiggle-config.yml",
    output:
        "results/Figure01/Figure01_prediction_ideal.slow5"
    conda:
        "../../envs/seq2squiggle.yml"
    log:
        "results/logs/run_seq2squiggle_Figure_01.log"
    threads: 8
    resources:
        disk_mb=50000,
        runtime=add_slack(100),
        mem_mb=35000,
        slurm="gpus=1",
        partition="zki", 
    shell:
        """
        resources/seq2squiggle/src/seq2squiggle/seq2squiggle.py predict --read-input --config {input.config} --model {input.model}  --ideal-event-length -1.0 --duration-sampling True --noise-sampling False --noise-std 0 --slow5 {output} {input.fasta}
        """


rule plot_Figure_01:
    input:
        events_tsv=rules.eventalign_valid_norm.output,
        slow5=rules.run_seq2squiggle_Figure_01.output,
        slow5_ideal=rules.run_seq2squiggle_Figure_01_ideal.output,
    output:
        "results/Figure01/Figure01.png"
    threads: 128
    conda:
        "../../envs/sam-stats.yml"
    log:
        "results/logs/sam-stats-plot.log",
    shell:
        """
        python workflow/scripts/Figure01_manuscript.py {input.events_tsv} {input.slow5} {input.slow5_ideal} {output}
        """

rule plot_sam_stats_groupplots:
    input:
        human_rm_seq2squiggle=rules.sam_stats_seq2squiggle.output,
        human_rm_squigulator=rules.sam_stats_squigulator.output,
        human_rm_exp=rules.sam_stats_experimental.output,
        human_gm_seq2squiggle=rules.sam_stats_seq2squiggle_gm.output,
        human_gm_squigulator=rules.sam_stats_squigulator_gm.output,
        human_gm_exp=rules.sam_stats_experimental_gm.output,
        fly_gm_seq2squiggle=rules.sam_stats_seq2squiggle_gm_fly.output,
        fly_gm_squigulator=rules.sam_stats_squigulator_gm_fly.output,
        fly_gm_exp=rules.sam_stats_experimental_gm_fly.output,
        ecoli_seq2squiggle=rules.sam_stats_seq2squiggle_gm_ecoli.output,
        ecoli_squigulator=rules.sam_stats_squigulator_gm_ecoli.output,
        ecoli_exp=rules.sam_stats_experimental_gm_ecoli.output,
    output:
        directory("results/group-plots/R10"),
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
        python workflow/scripts/basecalling-stats.py groupplot --outdir {output} \
        {input.human_rm_seq2squiggle} {input.human_rm_squigulator} {input.human_rm_exp} \
        {input.human_gm_seq2squiggle} {input.human_gm_squigulator} {input.human_gm_exp} \
        {input.fly_gm_seq2squiggle} {input.fly_gm_squigulator} {input.fly_gm_exp} \
        {input.ecoli_seq2squiggle} {input.ecoli_squigulator} {input.ecoli_exp} \
        --group_labels seq2squiggle,squigulator,experimental \
        --experiment_labels 'Human \nRead Mode,Human \nGenome Mode,D. melanogaster \nGenome Mode,E. coli \nGenome Mode'

        python workflow/scripts/basecalling-stats.py groupplotfigure --outdir {output} \
        {input.human_rm_seq2squiggle} {input.human_rm_squigulator} {input.human_rm_exp} \
        {input.human_gm_seq2squiggle} {input.human_gm_squigulator} {input.human_gm_exp} \
        {input.fly_gm_seq2squiggle} {input.fly_gm_squigulator} {input.fly_gm_exp} \
        {input.ecoli_seq2squiggle} {input.ecoli_squigulator} {input.ecoli_exp} \
        --group_labels seq2squiggle,squigulator,experimental \
        --experiment_labels 'Human \nRead Mode,Human \nGenome Mode,D. melanogaster \nGenome Mode,E. coli \nGenome Mode'
        """

rule plot_sam_stats_groupplots_R9:
    input:
        human_gm_seq2squiggle=rules.sam_stats_seq2squiggle_gm_R9.output,
        human_gm_squigulator=rules.sam_stats_squigulator_gm_R9.output,
        human_gm_exp=rules.sam_stats_experimental_gm_R9.output,
        human_gm_deepsimulator_CD=rules.sam_stats_deepsimulator_CI_gm_R9.output,
        human_gm_deepsimulator_CI=rules.sam_stats_deepsimulator_CD_gm_R9.output,
        sarscov_gm_seq2squiggle=rules.sam_stats_seq2squiggle_gm_sarscov_R9.output,
        sarscov_gm_squigulator=rules.sam_stats_squigulator_gm_sarscov_R9.output,
        sarscov_gm_exp=rules.sam_stats_experimental_gm_sarscov_R9.output,
        sarscov_gm_deepsimulator_CD=rules.sam_stats_deepsimulator_CD_gm_sarscov_R9.output,
        sarscov_gm_deepsimulator_CI=rules.sam_stats_deepsimulator_CI_gm_sarscov_R9.output,
    output:
        directory("results/group-plots/R9"),
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
        python workflow/scripts/basecalling-stats.py groupplotfigure --outdir {output} \
        {input.human_gm_seq2squiggle} {input.human_gm_squigulator} {input.human_gm_exp} \
        {input.human_gm_deepsimulator_CD} {input.human_gm_deepsimulator_CI} \
        {input.sarscov_gm_seq2squiggle} {input.sarscov_gm_squigulator} {input.sarscov_gm_exp} \
        {input.sarscov_gm_deepsimulator_CD} {input.sarscov_gm_deepsimulator_CI} \
        --group_labels 'seq2squiggle,squigulator,experimental,deepsimulator (context-dependent), deepsimulator (context-independent)' \
        --experiment_labels 'Human \nGenome Mode,SARS-CoV2 \nGenome Mode'
        """




rule publish_results_groupplots:
    input:
        groupplots_r10=rules.plot_sam_stats_groupplots.output,
        groupplots_r9=rules.plot_sam_stats_groupplots_R9.output,
        # figure01=rules.plot_Figure_01.output,
        config="config/config.yml",
        config_seq2squiggle="config/seq2squiggle-config.yml"
    output:
        directory("summary/group-plots/")
    log:
        "results/logs/publish_results_groupplots.log",
    shell:
        # Save the config file as well since we often run the different parts
        # of the benchmark at different times
        # cp -r {input.figure01} {output}
        """
        mkdir -p {output}
        cp -r {input.groupplots_r10} {output}
        cp -r {input.groupplots_r9} {output}
        cp {input.config} {output}
        cp {input.config_seq2squiggle} {output}
        """