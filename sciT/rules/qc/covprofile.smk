rule computeMatrix:
    input:
        bw = 'results/libraries/{library}/{name}.bw'
    output:
        density = 'results/libraries/{library}/{name}.density.tab.gz'
    threads: 12
    conda:
        "deeptools.yaml"
    params:
        gtf=config['reference']["annotation_gtf"]
    log:
        stdout="log/compute_matrix/{library}.{name}.stdout",
        stderr="log/compute_matrix/{library}.{name}.stderr"
    resources:
        mem_mb=lambda wildcards, attempt: 10000 * attempt, # Was benchmarked at 5041 megabytes
        runtime=lambda wildcards, attempt: f"{2030*attempt}s" # Was benchmarked at 685 seconds
    threads: 12
    shell:
        #"computeMatrix reference-point -S {input.bw} -R /references/ensembl/97/homo_sapiens/genes.bed --referencePoint TSS -a 3000 -b 3000 -out {output.density} -p {threads}"
        "computeMatrix scale-regions -S {input.bw} "
        "-R {params.gtf} "
        "-p {threads} "
        "--skipZeros --regionBodyLength 5000 -a 3000 -b 3000 "
        "-out {output.density} -p {threads} "
        "> {log.stdout} 2> {log.stderr}"

rule plotProfile:
    input:
        density = 'results/libraries/{library}/{name}.density.tab.gz'
    output:
        profile = 'results/libraries/{library}/{library}.{name}.plotProfile.tab',
        profile_plot = 'results/libraries/{library}/{library}.{name}.plotProfile.png'
    conda:
        "deeptools.yaml"
    log:
        stdout="log/plotProfile/{library}.{name}.stdout",
        stderr="log/plotProfile/{library}.{name}.stderr"
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt, # Was benchmarked at 3089 megabytes
        runtime=lambda wildcards, attempt: f"{600*attempt}s" # Was benchmarked at 44 seconds
    shell:
        "plotProfile -m {input.density} "
        "-o {output.profile_plot} "
        "--outFileNameData {output.profile} "
        "> {log.stdout} 2> {log.stderr}"
