rule samtools_coverage:
    input:
        bam = 'results/libraries/{library}/{prefix}.bam'
    output:
        report = 'results/libraries/{library}/{prefix}.csv'
    log:
        stderr="log/samtools_coverage/{library}_{prefix}.stderr"

    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt, # Was benchmarked at 90 megabytes
        runtime=lambda wildcards, attempt: f"{400*attempt}s" # Was benchmarked at 187 seconds
    shell:
        """samtools coverage {input.bam} | grep -v 'JH\\|GL\\|MU\\|KI' | awk 'NR<3{{print $0;next}}{{print $0| "sort --version-sort"}}' > {output.report} """

rule bam_index:
    input:
        bam = '{prefix}.bam'
    output:
        bam = '{prefix}.bam.bai'
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt, # Was benchmarked at 34 megabytes
        runtime=lambda wildcards, attempt: f"{300*attempt}s" # Was benchmarked at 22 seconds
    shell:
        "samtools index {input.bam} -@{threads}"
