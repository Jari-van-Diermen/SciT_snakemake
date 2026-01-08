
try:
    sciT_tagging_version = subprocess.check_output("sciT-tagger -v", shell=True)
except subprocess.CalledProcessError:
    raise ValueError('sciT-tagger seems not to be installed yet, check the manual section "Install DamID packages", and make sure the environment is activated')

rule sciT_tagging_se:
    input:
        #unsortedbam = "results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out.bam",
        unsortedbam = "results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out_htseq.bam",
        unsortedbam_index = "results/libraries/{library}/RNA_SE_Aligned.sortedByCoord.out_htseq.bam.bai",
        reference_fasta_path = config['reference']['fasta'],
        multiomes = "results/libraries/{library}/multiomes.yaml",
        meta = "results/libraries/{library}/meta.yaml"
    output:
        transcriptome = "results/libraries/{library}/transcriptome_se.bam",
        transcriptome_index = "results/libraries/{library}/transcriptome_se.bam.bai",
        statistics = "results/libraries/{library}/se_tagging.yaml"
        
    resources:
        mem_mb=lambda wildcards, attempt: 17409 * attempt ,  # Was benchmarked at 13391 megabytes
        runtime=lambda wildcards, attempt: f"{1512*attempt}s" # Was benchmarked at 1007 seconds
    params:
        gene_tag='XF', # use GN for using the STAR assigned counts, these are lower as they are stranded and stricter.
        sciT_tagging_version=sciT_tagging_version
    log:
        stdout="log/tag/{library}.stdout",
        stderr="log/tag/{library}.stderr"
    threads: 16
    shell:
        "sciT-tagger "
        "--threads {threads} "
        "-bam_in {input.unsortedbam} "
        "-reference_fasta_path {input.reference_fasta_path} "
        "-transcriptome_out {output.transcriptome} "
        "-statistics {output.statistics} "
        "-multiomes {input.multiomes} "
        "-gene_tag {params.gene_tag} "
        "-meta {input.meta} "
        "> {log.stdout} 2> {log.stderr}"

