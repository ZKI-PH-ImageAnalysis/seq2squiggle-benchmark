# Helpers to reduce duplication in the way we handle timeout, out of memory,
# and other failures
timeout_pre = """
        RETVAL=0;
        workflow/scripts/tout --dutmpdir -z {output.status} -k {params.max_runtime} -s {params.max_mem} -c -- bash -c '\
        """

# A "status" output file needs to be added for this to work properly
timeout_post = """
        echo "Return value: $RETVAL" >> {log}
        """

# Simplify adding a bit of slack to the SLURM resource limits, so that the tout
# script kills the job instead of SLURM
slack_multiplier = 1.1  # 10% slack


def add_slack(orig_limit):
    new_limit = round(orig_limit * slack_multiplier)
    return new_limit


def remove_slack_mem(wildcards, input, output, threads, resources):
    mem_limit = 1024 * round(resources.mem_mb / slack_multiplier)
    return mem_limit


def remove_slack_runtime(wildcards, input, output, threads, resources):
    runtime_limit = 60 * round(resources.runtime / slack_multiplier)
    return runtime_limit


rule save_config:
    input:
        "config/config.yml"
    output:
        "results/config.yml"
    shell:
        """
        cp {input} {output}
        """