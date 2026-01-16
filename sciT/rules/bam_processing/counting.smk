
# Associate reads to genes
rule transcriptome_htseq:
    input:
        ref=config['reference']['fasta'],
        bam="results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out.bam",
        annotations=config['reference']["annotation_gtf"]

    output:
        #counts = temp("results/libraries/{library}/transcriptome_se.loom")
        bam=temp("results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out_htseq.bam")
    log:
        stdout="log/se_counting/{library}.stdout",
        stderr="log/se_counting/{library}.stderr"
    conda:
        "htseq.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,  # Was benchmarked at 98 megabytes
        runtime=lambda wildcards, attempt: f"{12926*attempt}s"  # Was benchmarked at 8616 seconds
    shell:
        "htseq-count {input.bam} "
        "{input.annotations} "
        "--order pos "
        "-s no " # the data is not stranded
        "-i gene_id "
        "-t gene "
        "-o {output.bam} "
        "-p BAM "
        "> {log.stdout} 2> {log.stderr}"




rule sciT_cell_filter:
    input:
        transcriptome_bam = "results/libraries/{library}/transcriptome_se.bam"
    params:
        cell_min_transcriptome_count = (config.get('sciT', {}).get('cell_min_transcriptome_count', 0))
    output:
        transcriptome_output_bam = "results/libraries/{library}/transcriptome_se.cell_filtered.bam",
        rg_filetx = temp("results/libraries/{library}/transcriptome_se.cell_filtered.rg.txt"),
        cell_listtx = "results/libraries/{library}/transcriptome_se.cell_filtered.samples.txt"
    log:
        stdout="log/sciT_cell_filter/{library}.stdout",
        stderr="log/sciT_cell_filter/{library}.stderr"
    shell:
        "sciT-cell-filter "
        "--cell_min_transcriptome_count {params.cell_min_transcriptome_count} "
        "--transcriptome_input_bam {input.transcriptome_bam} "
        "--transcriptome_output_bam {output.transcriptome_output_bam} "
        "> {log.stdout} 2> {log.stderr}"




if 'hybrid' in config['reference']:
    # Create count files when data  is allele specific
    rule tx_bam2loom:
        input:
            bam= "results/libraries/{library}/transcriptome_se.cell_filtered.bam",
            annotations = config['reference']['annotation_gtf']
        output:
            loom = "results/libraries/{library}/transcriptome_se.loom"
        log:
            stdout="log/epibam2loom/{library}_transcriptome_se.stdout",
            stderr="log/epibam2loom/{library}_transcriptome_se.stderr"
        threads:22
        params:
            gene_tag='XF',
            version = subprocess.check_output("transcriptomebam-to-loom --version", shell=True)
        resources:
            mem_mb=lambda wildcards, attempt: 15000 * attempt,
            runtime=lambda wildcards, attempt: f"{7200*attempt}s"  
        shell:
            "transcriptomebam-to-loom "
            "{input.bam} "
            "--threads {threads} "
            "--gene_tag {params.gene_tag} "
            "--loom_out {output.loom} "
            "--annotations {input.annotations} "
            "--allele_specific "
            "> {log.stdout} 2> {log.stderr}"

else:
    rule tx_bam2loom:
        input:
            bam="results/libraries/{library}/transcriptome_se.cell_filtered.bam",
            annotations = config['reference']['annotation_gtf']
        output:
            loom = "results/libraries/{library}/transcriptome_se.loom"
        log:
            stdout="log/epibam2loom/{library}_transcriptome_se.stdout",
            stderr="log/epibam2loom/{library}_transcriptome_se.stderr"
        threads:22
        params:
            gene_tag='XF',
            version = subprocess.check_output("transcriptomebam-to-loom --version", shell=True)
        resources:
            mem_mb=lambda wildcards, attempt: 15000 * attempt,  # Was benchmarked at 6062 megabytes
            runtime=lambda wildcards, attempt: f"{7200*attempt}s"  # Was benchmarked at 194 seconds
        shell:
            "transcriptomebam-to-loom "
            "{input.bam} "
            "--threads {threads} "
            "--gene_tag {params.gene_tag} "
            "--loom_out {output.loom} "
            "--annotations {input.annotations} "
            "> {log.stdout} 2> {log.stderr}"





