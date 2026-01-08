rule fastq_screen:
    input:
        fastq = "results/libraries/{library}/trimmed.{suffix}.fastq.gz",
        config = "fastq_screen.conf"
    output:
        "results/libraries/{library}/trimmed.{suffix}_screen.txt"
    log:
        stdout = "log/fastq_screen/{library}_{suffix}.stdout",
        stderr = "log/fastq_screen/{library}_{suffix}.stderr"

    conda:
        "fastq_screen.yaml"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt, # Was benchmarked at 21 megabytes, which is surely wrong
        runtime=lambda wildcards, attempt: f"{1085*attempt}s" # Was benchmarked at 722 seconds
    shell:
        "fastq_screen --aligner bowtie2 "
        "--threads {threads} "
        "--outdir results/libraries/{wildcards.library} "
        "--conf {input.config} "
        "--force {input.fastq} "
        "> {log.stdout} 2> {log.stderr}"

