rule index_STAR:
    input:
        ref=config['reference']['fasta']
    output:
        chlen = "star_genome_index/chrLength.txt"
    conda:
        "star.yaml"
    log:
        stdout="log/STAR_INDEX.stdout",
        stderr="log/STAR_INDEX.stderr"
    threads: 128

    resources:
        mem_mb=lambda wildcards, attempt: 71804 * attempt ,  # Was benchmarked at 27616 megabytes
        runtime=lambda wildcards, attempt: f"{3000*attempt}s" # Was benchmarked at 703 seconds
    shell:
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir star_genome_index "
        "--genomeFastaFiles {input.ref} "
        "> {log.stdout} 2> {log.stderr}"

if 'hybrid' not in config['reference']:
    rule STAR_PE:
        input:
            starIndex = "star_genome_index/chrLength.txt",
            r1 = "results/libraries/{library}/trimmed.R1.fastq.gz",
            r2 = "results/libraries/{library}/trimmed.R2.fastq.gz",
        output:
            bam = "results/libraries/{library}/RNA_PE_Aligned.sortedByCoord.out.bam",
            counts = "results/libraries/{library}/RNA_PE_ReadsPerGene.out.tab",
            logfile = "results/libraries/{library}/RNA_PE_Log.final.out",
            pass1_log = temp("results/libraries/{library}/RNA_PE__STARpass1/Log.final.out"),
            pass1_sj = temp("results/libraries/{library}/RNA_PE__STARpass1/SJ.out.tab"),
            pass1folder = temp(directory("results/libraries/{library}/RNA_PE__STARpass1")),
            g1folder = temp(directory("results/libraries/{library}/RNA_PE__STARgenome"))

        threads: 32
        log:
            stdout="log/STAR_PE/{library}.stdout",
            stderr="log/STAR_PE/{library}.stderr"
        params:
            gtf=config['reference']["annotation_gtf"]

        conda:
            "star.yaml"
        shell:
            "STAR "
            "--runThreadN {threads} "
            "--runMode alignReads "
            "--outSAMunmapped Within "
            "--genomeDir star_genome_index "
            "--readFilesIn {input.r1} {input.r2} "
            "--outFileNamePrefix results/libraries/{wildcards.library}/RNA_PE_ "
            "--sjdbGTFfile {params.gtf} "
            "--quantMode GeneCounts "
            "--sjdbGTFtagExonParentGene gene_name "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI nM AS GX GN jM jI "
            "--readFilesCommand zcat "
            "--twopassMode Basic"
            "> {log.stdout} 2> {log.stderr}"
else:
    rule STAR_PE_WASP:
        input:
            starIndex = "star_genome_index/chrLength.txt",
            r1 = "results/libraries/{library}/trimmed.R1.fastq.gz",
            r2 = "results/libraries/{library}/trimmed.R2.fastq.gz",
            varVCFfile = "references/hybrid/snps.vcf"
        output:
            bam = "results/libraries/{library}/RNA_PE_Aligned.sortedByCoord.out.bam",
            counts = "results/libraries/{library}/RNA_PE_ReadsPerGene.out.tab",
            logfile = "results/libraries/{library}/RNA_PE_Log.final.out",
            pass1_log = temp("results/libraries/{library}/RNA_PE__STARpass1/Log.final.out"),
            pass1_sj = temp("results/libraries/{library}/RNA_PE__STARpass1/SJ.out.tab"),
            pass1folder = temp(directory("results/libraries/{library}/RNA_PE__STARpass1")),
            g1folder = temp(directory("results/libraries/{library}/RNA_PE__STARgenome"))

        threads: 32
        log:
            stdout="log/STAR_PE/{library}.stdout",
            stderr="log/STAR_PE/{library}.stderr"
        params:
            gtf=config['reference']["annotation_gtf"]

        conda:
            "star.yaml"
        shell:
            "STAR "
            "--runThreadN {threads} "
            "--runMode alignReads "
            "--outSAMunmapped Within "
            "--genomeDir star_genome_index "
            "--readFilesIn {input.r1} {input.r2} "
            "--outFileNamePrefix results/libraries/{wildcards.library}/RNA_PE_ "
            "--sjdbGTFfile {params.gtf} "
            "--quantMode GeneCounts "
            "--sjdbGTFtagExonParentGene gene_name "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI nM AS GX GN rB vA vG vW jM jI "
            "--readFilesCommand zcat "
            "--varVCFfile {input.varVCFfile} "
            "--twopassMode Basic"
            "> {log.stdout} 2> {log.stderr}"


if 'hybrid' not in config['reference']:
    # Use single end mapping of the second mate for "normal" CS2 barcoded reads, the first read is pretty much useless
    rule STAR_transcriptome_single_end:
        input:
            starIndex='star_genome_index/chrLength.txt',
            r2 = "results/libraries/{library}/RNA_barcoded.se.fastq.gz"

        output:
            bam="results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out.bam",
            counts="results/libraries/{library}/RNA_SE_ReadsPerGene.out.tab",
            logfile="results/libraries/{library}/RNA_SE_Log.final.out" #,
            #pass1_log = temp("results/libraries/{library}/RNA_SE__STARpass1/Log.final.out"),
            #pass1_sj = temp("results/libraries/{library}/RNA_SE__STARpass1/SJ.out.tab"),
            #pass1folder = temp(directory("results/libraries/{library}/RNA_SE__STARpass1")),
            #g1folder = temp(directory("results/libraries/{library}/RNA_SE__STARgenome"))
        threads: 32
        log:
            stdout="log/STAR_SE/{library}.stdout",
            stderr="log/STAR_SE/{library}.stderr"
        params:
            gtf=config['reference']["annotation_gtf"]

        conda:
            "star.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 80000 * attempt ,  # Was benchmarked at 31708 megabytes
            runtime=lambda wildcards, attempt: f"{3000*attempt}s" # Was benchmarked at 534 seconds
        shell:
            "STAR "
            "--runThreadN {threads} "
            "--runMode alignReads "
            "--genomeDir star_genome_index "
            "--readFilesIn {input.r2} "
            "--outSAMunmapped Within "
            "--outFileNamePrefix results/libraries/{wildcards.library}/RNA_SE_ "
            "--sjdbGTFfile {params.gtf} "
            "--quantMode GeneCounts "
            "--sjdbOverhang 74 "
            "--sjdbGTFtagExonParentGene gene_id "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI nM AS GX GN MD jM jI "
            "--readFilesCommand zcat " #"--twopassMode Basic"
            "> {log.stdout} 2> {log.stderr}"
else:
    rule STAR_transcriptome_single_end_WASP:
        input:
            starIndex='star_genome_index/chrLength.txt',
            r2 = "results/libraries/{library}/RNA_barcoded.se.fastq.gz",
            varVCFfile = "references/hybrid/snps.vcf"

        output:
            bam="results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out.bam",
            counts="results/libraries/{library}/RNA_SE_ReadsPerGene.out.tab",
            logfile="results/libraries/{library}/RNA_SE_Log.final.out" #,
            #pass1_log = temp("results/libraries/{library}/RNA_SE__STARpass1/Log.final.out"),
            #pass1_sj = temp("results/libraries/{library}/RNA_SE__STARpass1/SJ.out.tab"),
            #pass1folder = temp(directory("results/libraries/{library}/RNA_SE__STARpass1")),
            #g1folder = temp(directory("results/libraries/{library}/RNA_SE__STARgenome"))
        threads: 32
        log:
            stdout="log/STAR_SE/{library}.stdout",
            stderr="log/STAR_SE/{library}.stderr"
        params:
            gtf=config['reference']["annotation_gtf"]

        conda:
            "star.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: 80000 * attempt ,  # Was benchmarked at 33813 megabytes
            runtime=lambda wildcards, attempt: f"{3000*attempt}s" # Was benchmarked at 658 seconds
        shell:
            "STAR "
            "--runThreadN {threads} "
            "--runMode alignReads "
            "--genomeDir star_genome_index "
            "--readFilesIn {input.r2} "
            "--outSAMunmapped Within "
            "--waspOutputMode SAMtag "
            "--outFileNamePrefix results/libraries/{wildcards.library}/RNA_SE_ "
            "--sjdbGTFfile {params.gtf} "
            "--quantMode GeneCounts "
            "--sjdbOverhang 74 "
            "--sjdbGTFtagExonParentGene gene_id "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI nM AS GX GN MD rB vA vG vW jM jI "
            "--readFilesCommand zcat " #"--twopassMode Basic"
            "--varVCFfile {input.varVCFfile} "
            "> {log.stdout} 2> {log.stderr}"



