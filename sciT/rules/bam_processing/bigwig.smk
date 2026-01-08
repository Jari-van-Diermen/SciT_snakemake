rule bigwig:
    input:
        "results/libraries/{library}/{name}.bam",
        "results/libraries/{library}/{name}.bam.bai"
    output:
        "results/libraries/{library}/{name}.bw"
    conda:
        "deeptools.yaml"

    params:
        binsize=50

    log:
        stdout="log/bigwig/{library}_{name}.stdout",
        stderr="log/bigwig/{library}_{name}.stderr"

    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt, # Was benchmarked at 950 megabytes
        runtime=lambda wildcards, attempt: f"{300*attempt}s" # Was benchmarked at 159 seconds
    shell:
        "bamCoverage "
        "-b {input[0]} "
        "-o {output} "
        "-p {threads} "
        "--ignoreDuplicates --samFlagExclude 1536 "
        "--normalizeUsing CPM "
        "--binSize {params.binsize} > {log.stdout} 2> {log.stderr}"
