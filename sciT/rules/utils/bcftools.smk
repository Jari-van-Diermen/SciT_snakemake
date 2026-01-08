if 'hybrid' in config['reference']:
    # Generate a VCF file for haplotype A and B
    rule bcftools_subset_sample:
        input:
            vcf_file = config['reference']['hybrid']['multisample_vcf_file']
        output:
            vcf_file = "references/{haplotype}/snps.vcf.gz",  # Haplotype can be alleleA or alleleB
            label_file = "references/{haplotype}/README.txt",
            vcf_index = "references/{haplotype}/snps.vcf.gz.csi"
        log:
            stdout = "log/bcftools_subset_sample/{haplotype}.stdout",
            stderr = "log/bcftools_subset_sample/{haplotype}.stderr"
        conda:
            "bcftools.yaml"
        threads: 3
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * attempt,
            runtime=lambda wildcards, attempt: f"{2082*attempt}s"
        params:
            sample_name= lambda wildcards: config['reference']['hybrid'][wildcards.haplotype]
        wildcard_constraints:
            haplotype='alleleA|alleleB'
        shell:
            "echo 'This is data from haplotype {params.sample_name}' > {output.label_file} && bcftools view -s {params.sample_name} {input.vcf_file} --threads {threads} -O z -o {output.vcf_file} --write-index"
            " > {log.stdout} 2> {log.stderr}"
        

    # Generate a consensus fasta file for haplotype A and B
    rule bcftools_consensus:
        input:
            vcf_file = "references/{haplotype}/snps.vcf.gz",
            fasta_file = config['reference']['fasta']
        output:
            fasta_file = "references/{haplotype}/consensus.fasta"
        conda:
            "bcftools.yaml"
        log:
            stderr = "log/bcftools_consensus/{haplotype}.stderr"
        wildcard_constraints:
            haplotype='alleleA|alleleB'
        resources:
            mem_mb=lambda wildcards, attempt: 1024 * attempt,
            runtime=lambda wildcards, attempt: f"{1000*attempt}s"
        shell:
            "cat {input.fasta_file} | bcftools consensus {input.vcf_file} > {output.fasta_file} "
            "2> {log.stderr}"
