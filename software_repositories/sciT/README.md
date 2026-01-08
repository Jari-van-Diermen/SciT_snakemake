# Tools for sciT analysis

This is a set of tools for single-cell transcriptome analysis using single-cell combinatorial indexing


## sciT-tagger
Deduplicates transcriptome data based on the UMI, which is obtained using the `RX` sam tag.
Make sure the input file is sorted by query name! The outputs are sorted by coordinate.
When library-meta information is supplied, this script encodes library-meta information in DS field of the read groups of the output bam files.

```
usage: sciT-tagger [-h] -bam_in BAM_IN [-transcriptome_out TRANSCRIPTOME_OUT] -reference_fasta_path REFERENCE_FASTA_PATH
                    [-statistics STATISTICS] [-multiomes MULTIOMES] [-meta META] [--threads THREADS] [-gene_tag GENE_TAG]


options:
  -h, --help            show this help message and exit
  -bam_in BAM_IN        Path to input transcriptome bam file (default: None)
  -transcriptome_out TRANSCRIPTOME_OUT
                        Path to output transcriptome bam file
  -reference_fasta_path REFERENCE_FASTA_PATH
                        Path to reference fasta file (default: None)
  -statistics STATISTICS
                        Path to statistics yaml (default: None)
  -multiomes MULTIOMES  Multiomes file (yaml) (default: None)
  -meta META            Meta file (yaml), should contain a library-meta key, containing key:value pairs for the metadata of the different libraries.
  --threads THREADS     Multithread by splitting by chromosome, this can only be done when the input is single end and coordinate sorted (default: 1)
  -gene_tag GENE_TAG    Tag present in the input bam file which associates reads to a gene, for example GN for STAR, and XF for HTSEQ
                        (default is GN, however, the snakemake rule sciT_tagging_se sets it to XF).
  ```

## transcriptomebam-to-loom
Create a transcriptome loom file from a BAM file generated using sciT-tagger.
```
usage: transcriptomebam-to-loom [-h] [--loom_out LOOM_OUT] [--gene_tag GENE_TAG] [--annotations ANNOTATIONS] [--threads THREADS] bam_in

Extract an loom file from a tagged BAM file

positional arguments:
  bam_in                Path to input bam files

options:
  -h, --help            show this help message and exit
  --loom_out LOOM_OUT   Write loom file here (default: None)
  --gene_tag GENE_TAG   Tag which determines the gene (default: GN)
  --annotations ANNOTATIONS
                        Annotation file, when supplied this is used to also encode the GENE NAME in the loom file (default: None)
  --threads THREADS     Worker threads (default: None)

```

## sciT-qc

Run quality control and clustering on supplied loom files, using provided thresholds and meta information

```
usage: sciT-qc [-h] [--tx TX] [--tx_count_threshold TX_COUNT_THRESHOLD] [--config_file CONFIG_FILE]
                [--tx_clustering_output_path TX_CLUSTERING_OUTPUT_PATH] [--qc_csv QC_CSV]


options:
  -h, --help            show this help message and exit

Inputs:
  --tx TX               Path to transcriptome loom file, multiple files can be supplied, the selected CONDITION COLUMNS need to be encoded in this file (default:
                        None)

Thresholds:
  --tx_count_threshold TX_COUNT_THRESHOLD
                        Minimum required number of unique RNA reads to pass QC (default: 5000)

Conditions:
  --config_file CONFIG_FILE
                        Path to config.yaml file which has a key conditions: attributes (default: None)

Outputs:
  --tx_clustering_output_path TX_CLUSTERING_OUTPUT_PATH
                        path to transcriptome clustering plot output (default: None)
  --qc_csv QC_CSV       path to qc csv file (default: None)

```
