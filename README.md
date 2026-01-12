# snakemake pipeline for sciT analysis

This repository contains the scripts and instructions you will need to perform a sciT analysis.

This snakemake workflow represents an adaptation of an existing in-house snakemake workflow, which allows
it to preprocess single-cell transcriptomics data that was generated using a single-cell combinatorial
indexing approach (sciT).

Below you will find instructions for the [initial setup](#initial-setup) of the sciT working environment, and also the
required steps you perform for [each sciT workflow](#execution-of-scit-workflow).

## initial setup

### Install miniforge3

The first time you run the sciT-snakemake pipeline, you have to prepare an environment that can prepare and run the neccesary code.

You can use conda for this. I recommend you to use *Miniforge3* as a conda package manager. It can be installed like this:

```bash
# Use wget to download the installer
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# Or use curl to download the installer
#curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

# Run the installer
bash Miniforge3-$(uname)-$(uname -m).sh
```

> [!NOTE]
> Using the installer like above should run it in interactive mode, allowing you to indicate the install directory for Miniforge3.

> [!NOTE]
> The installer will prompt you to initialize conda with your shell. Select yes if you want conda to activate automatically.

> [!IMPORTANT]
> If installing miniforge3 in an HPC, do **NOT** install it in your home directory. It will quickly run out of storage space for your environments.
>

### Create base conda environment

Now you can create the base conda environment that will be used for running the sciT workflow.

After conda is activated, run the following code to create the base environment.

```bash
conda create -c conda-forge -c bioconda -n sciT_snakemake -y snakemake samtools p7zip gcc cxx-compiler zlib openjdk snakemake-executor-plugin-slurm mamba wget git
```

### Check reference repositories

For the sciT workflow to work, you need proper reference files. These can be placed in a `reference_repositories` directory (When executing the sciT workflow we will link to the exact location of the `reference_repositories` directory, so it can be anywhere you like). This should include:

- **contaminants:** A directory with *bowtie2* genome indexes and fasta files for potential contaminants to test for using *fastq-screen*. for example:

  ```txt
  |____reference_repositiories
  | |____contaminants
  | | |____cutibacterium.1.bt2
  | | |____cutibacterium.2.bt2
  | | |____cutibacterium.3.bt2
  | | |____cutibacterium.4.bt2
  | | |____cutibacterium.rev.2.bt2
  | | |____cutibacterium.rev.1.bt2
  | | |____cutibacterium.fasta
  ```

- Other directories containing genome files and gtf annotations for species you want to map your read onto. For example, here is a directory for **mm39 references**:

  ```txt
  |____reference_repositiories
  | |____mm39
  | | |____annotations.gtf
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.1.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.2.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.3.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.4.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.rev.1.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.rev.2.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.fai
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.dict
  | | |____mgp_REL2021_snps.vcf.gz.csi
  | | |____mgp_REL2021_snps.vcf.gz.tbi
  | | |____mgp_REL2021_snps.vcf.gz
  ```

  Here, there are

  - **bowtie index files:** for including *mm39* in the fastq-screen report.
  - **fasta file:** for building the index for the [STAR aligner](https://pubmed.ncbi.nlm.nih.gov/23104886/).
  - **annotations.gtf:** For gene counting by [HTseq](https://htseq.readthedocs.io/en/latest/).

> [!NOTE]
> For the references to be included in the *fastq-screen* report, they need to be included in the `fastq_screen.conf` file.
> 
>     DATABASE	Mouse	reference_repositories/mm39/Mus_musculus.GRCm39.dna_sm.primary_assembly
>     DATABASE	Cutibacterium_Acnes	reference_repositories/contaminants/cutibacterium
  
When you run the sciT worklfow, you will select *one reference fasta file* and one *GTF annotation file*.
This will be selected in the *config.yaml* file (more on that later).

### Install sciT workflow packages

In order to make the sciT workflow work, we need to install some additional packages. These are located in the `software_repositories` directory in this github repository. To install these packages in the **sciT_snakemake** conda environment, first activate the environment:

```bash
conda activate sciT_snakemake
```

Clone this github repository

```bash
# clone repository
git clone git@github.com:Suirotras/SciT_snakemake.git
# move into repository
cd SciT_snakemake
```

Then install the packages into the conda environment using **pip**.

```bash
pip install ./software_repositories/* --upgrade
```

This installation only needs to be performed once, after which this **sciT_snakemake** environment can be used for different sciT workflows.

You can now delete the cloned github repository, as you only used it here to install the sciT workflow packages.

## Execution of sciT workflow

### Activate the base environment.

```bash
conda activate sciT_snakemake
```

### Clone repository

Clone this github repository to the location where you want to preprocess the sciT fastq.gz files.

```bash
cd path/to/desired/location
# clone repository
git clone git@github.com:Suirotras/SciT_snakemake.git
# move into repository
cd SciT_snakemake
cd sciT
```

> [!IMPORTANT]
> If this is the FIRST time you are running this workflow, you NEED to first install the packages in the `software_repositories` directory.
> Check out the [**Install sciT workflow packages**](#install-scit-workflow-packages) section first before continuing.
>
> If this was not the first time running this workflow, you can safely remove the `software_repositories` directory, as it won't be needed.

### Create a symlink to the reference repository directory

You now need to tell where the `reference_repositories` live, and subsequently create a symbolic link to that directory. If you do not do this, then the sciT workflow will not be able to find the references.

```bash
ln -s ../reference_repositories reference_repositories
```

### Modify config.yaml

You will need to configure the configuration file that tells the sciT workflow which options to use and which files to process. Here is an example of a correct `config.yaml`:

```yaml
reference:
  # Human Hg38:
  # fasta: reference_repositories/HG38.p14/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
  # annotation_gtf: reference_repositories/HG38.p14/annotations.gtf
  # For mouse mm10 use:
  fasta: reference_repositories/mm10_pool_lab/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC.fa
  annotation_gtf: reference_repositories/mm10_pool_lab/annotations.gtf
  # hybrid: # Comment or remove this section to not perform allele specific analysis
  #   multisample_vcf_file: reference_repositories/mm39/mgp_REL2021_snps.vcf.gz
  #   alleleA: 129S1_SvImJ # Corresponds to a sample name in the VCF file above
  #   alleleB: CAST_EiJ

QC: # The quality control thresholds used for the quality control report
  transcriptome_count_threshold: 0
  transcriptome_genes_threshold: 1000
  report_top_samples: 100 # Number of best cells/samples to include in the report, based on the number of unique sciT reads
  report_bottom_samples: 100 # Number of worst cells to include in the report, based on the number of unique sciT reads

sciT:
  cell_min_transcriptome_count: 25 # Minimum number of transcriptome molecules for a cell to be counted

conditions:
  attributes: # Attributes used to determine the condition (label) of each sample
  - number
  - library
  colors: # Predetermined colors for attribute values:
    number:
      1K: 0000ff
      10K: abc8a4
    library:
      i31: ffcc00

# Add library data here by creating a libraries section, see documentation for more details
libraries:
  i31-sample1-1K:
    type: sciT
    demultiplexer: sciT
    transcriptome-paired-end-fastq-files:
      R1:
        - i31-sample1_L001_R1_001.fastq.gz
        - i31-sample1_L002_R1_001.fastq.gz
      R2:
        - i31-sample1_L001_R2_001.fastq.gz
        - i31-sample1_L002_R2_001.fastq.gz
    library-meta:
      number: 1K
      library: i31
```

`config.yaml` consists of the following fields:

- **reference**

  - **fasta:** genome file from which STAR index will be created. Reads will be mapped to this STAR index.
  - **annotation_gtf:** GTF file that will be used to count transcriptomic features.
  - **hybrid:** Uncomment to perform an allele-specific analysis
    - **multisample_vcf_file:** vcf file needed for allele-specific analysis.
                                Uncomment to perform an allele-specific analysis.
    - **alleleA:** sample name in vcf file of allele A. Uncomment to perform an allele-specific analysis.
    - **alleleB:** sample name in vcf file of allele B. Uncomment to perform an allele-specific analysis.

- **QC:** Settings that will influence the QC report

- **sciT**
  - **cell_min_transcriptome_count:** A cell filter. Will filter out cells with less counts than the number indicated here.

- **conditions**
  - **attributes:** These attributes, in conjunction with the values in the *libraries/library-meta* fields, are used to select colors that will be used in
                    the quality report.
    - For example, a **seaborn clustermap** that is generated will annotate the identified cells based on the metadata in its *library-meta* field.
      In the case of the *config.yaml* above, the clustermap will have two annotation colorbars, one for *library* and one for *number* (i.e. number of
      cells used to generate the library).  
      - **Library colorbar**: cells from library i31 are *Tangerine yellow*
        (i.e. `#ffcc00` <span style="display:inline-block;width:12px;height:12px;background-color:#ffcc00;border:1px solid #000;"></span>).
      - **number colorbar**: cells from a library prepared with 1K cells are *Blue*
        (i.e. `#0000ff` <span style="display:inline-block;width:12px;height:12px;background-color:#0000ff;border:1px solid #000;"></span>),
        while cells from a library prepared with 10K cells are *Green*
        (i.e. `#abc8a4` <span style="display:inline-block;width:12px;height:12px;background-color:#abc8a4;border:1px solid #000;"></span>).

        Resulting in a plot seen below, as one library (i31) with 1K cells was processed:

        <img src="./figures/transcriptome_clustering_mqc.png" alt="example_clustermap" width="500">


- **libraries:** Here you put the information about individual libraries that will be processed.

  - **i31-sample1-1K:** Name of your library identifier.
    - **type:** NEEDS to be **sciT** for every library entry.
    - **demultiplexer:** NEEDS to be **sciT** for every library entry.
    - **transcriptome-paired-end-fastq-files:** Indicate the paired end fastq-files that represent this library and should be processed.
    - **library-meta:** This indicates the library-wide metadata that defines this library. Reads from these libraries will be annotated with this meta-data.

      - *Bam files* for this library will have read groups annotated with this metadata:

        `@RG	ID:... ... DS:BC:...;number:1K;library:i31`

      - *Loom files* for this library will have this *library-meta* data in its *Cell metadata* field.

### Modify fastq_screen.conf

In order for **fastq-screen** to work correctly, make sure the repositories are correctly listed in the `fastq_screen.conf`. Check this for both the
contaminants and references you are using for read mapping. This allows **fastq-screen** to find the correct bowtie genome indexes.

### Execute the sciT workflow

With all this set up, you can start the sciT workflow in order to process the sciT data.

#### for usage locally (i.e. your own laptop)

Adjust the core count and memory size

```bash
snakemake --cores 32 --resources mem_mb=62000 -k -p --rerun-incomplete --nt --use-conda --conda-prefix environments
```

#### for usage on SLURM cluster

The `--cores 60` and `mem_mb=100000` might not be the most ideal settings, You can probably change it to lower values.

```bash
snakemake --executor slurm --cores 60 --resources mem_mb=100000 -k -p --rerun-incomplete --jobs 20 --restart-times 3 --use-conda --conda-prefix environments
```

> [!NOTE]
> Use the `screen` utility to run the sciT workflow. This makes sure the sciT workflow can keep running when you to disconnect from the HPC.


## sciT-snakemake expected output

After the *sciT_snakemake* workflow, the output will be in a new subdirectory named `results`. This directory contains the following:

- **libraries:** per-library output files. Some notable output files are:
  - *`RNA_barcoded.se.fastq.gz`*: fastq file after barcode identification and read trimming.
  - *`RNA_SE_Aligned.sortedByCoord.out.bam`*: Bam created by mapping `RNA_barcoded.se.fastq.gz` with STAR.
  - *`transcriptome_se.bam`*: Bam file created by using `RNA_SE_Aligned.sortedByCoord.out.bam` to assign genes to reads.
                              Furthermore, reads are deduplicated here.
  - *`transcriptome_se.cell_filtered.bam"`*: Bam file after cell-filtering using the **sciT/cell_min_transcriptome_count**
                                            field in the *config.yaml*. Simply put, all barcodes with fewer UMI-counts than the
                                            value indicated in the **sciT/cell_min_transcriptome_count** field are excluded here.
  - *`transcriptome_se.loom`*: Loom file created by converting `transcriptome_se.cell_filtered.bam` to loom format.
- **QC:** quality-control data:
  - *`qc_status_per_sample.csv`*: comma-separated file indicating whether a barcode passes the QC threshold (e.g. the minimum
                                  UMI-counts indicated in the **sciT/cell_min_transcriptome_count** in the *config.yaml*).
  - *`transcriptome_clustering_mqc.png`*: Clustermap of the processed libraries.
- **report:**
  - *`qc_report.html`*: The MultiQC report.

## Using sciT-snakemake output

For further analyses of the processed libraries by the *sciT-snakemake* workflow, the `transcriptome_se.loom` files are the most convenient to use.

Here, I provide a convenient way of loading sciT-snakemake-generated loom files into R.

### Loading sciT-snakemake data in R

I provide a companion R package for this sciT-snakemake workflow. This package is named *LoadSciLooms*, and can be found on
github [here](https://github.com/Suirotras/LoadSciLooms). The package has convenient functions to load loom files as Seurat objects.

To use *LoadSciLooms*, simply install the package using R devtools.

```R
library(devtools)
devtools::install_github("Suirotras/LoadSciLooms", ref = "master")
```

After installing *LoadSciLooms*, the `LoomAsSeurat` function can be used to load a loom file as a seuratobject.

the *sciT-snakemake* workflow creates loom files with different metadata variables:

- **name:** metadata variable referring to gene symbols (e.g. Il12rb1)
- **Gene:** metadata variable referring to ENSEMBL gene IDs (e.g. ENSMUSG00000000791)

Thus, we can create a seuratobject using the gene symbols:

```R
# Load as seuratobject
Loom_path <- "path/to/libraries/i31-sample1-1K/transcriptome_se.loom"
sciT_seurat <- LoomAsSeurat(Loom_path, matrix_rowname_col = "name",
                            resolve_duplicates = TRUE,
                            gmm_cell_calling = FALSE)$seurat
```

or by using the ENSEMBL gene IDs:

```R
# Load as seuratobject
Loom_path <- "path/to/libraries/i31-sample1-1K/transcriptome_se.loom"
sciT_seurat <- LoomAsSeurat(Loom_path, matrix_rowname_col = "Gene",
                            resolve_duplicates = TRUE,
                            gmm_cell_calling = FALSE)$seurat
```

The *LoomAsSeurat* function also has extra functionalities:

- Determining which barcodes represent 'real' cells using gaussian mixed modeling (gmm).
  The `gmm_call_calling` parameter activates this behaviour.
- Automatically resolves duplicate feature names, which is neccessary when you create the
  Seurat object with gene symbols. This is controlled by the `resolve_duplicates` parameter.

Furthermore, multiple loom files can be opened as a single merged Seurat object using the
`MultiLoomAsSeurat` function:

For further information, check out the documentation of the *LoadSciLooms*.

## Issues

### Out of memory

It is possible that for some rules, not enough memory is allocated for your specific job.

For example, an `OUT_OF_MEMORY` error occured for me in the rule that concatenates fastq files.
In the snakemake log file, that looks like this:

```bash
Error in rule sciT-RNA-concat-i31-sample1-1K-R2:
    message: SLURM-job '45153759' failed, SLURM status is: 'OUT_OF_MEMORY'. For further error details see the cluster/cloud log and the log files of the involved rule(s).
    jobid: 13
    input: /path/to/subsampled/i31-sample1-sciT_L001_R2_001.fastq.gz, /path/to/subsampled/i31-sample1-sciT_L002_R2_001.fastq.gz
    output: results/libraries/i31-sample1-1K/raw.RNA.R2.fastq.gz
    log: log/concat/i31-sample1-1K_R2.stderr, /path/to/.snakemake/slurm_logs/rule_sciT-RNA-concat-i31-sample1-1K-R2/45153759.log (check log file(s) for error details)
    shell:
        cat /path/to/subsampled/i31-sample1-sciT_L001_R2_001.fastq.gz /path/to/subsampled/i31-sample1-sciT_L002_R2_001.fastq.gz > results/libraries/i31-sample1-1K/raw.RNA.R2.fastq.gz 2> log/concat/i31-sample1-1K_R2.stderr
        (command exited with non-zero exit code)
    external_jobid: 45153759
```

In these cases, you should find the problematic rule, which should be located in the `snakefile` or one of the
`rules/*/*.smk` files. This problematic rule was located in the `snakefile`.

After finding the problematic rule, find its **resources** directive, which looks something like this:

```bash
resources:
  mem_mb=lambda wildcards, attempt: attempt * 100
```

This essentially means that the memory (i.e. `mem_mb`) allocated to this rule on its first attempt is 100MB.
On the second attempt, the allocated memory will be 200MB.

The `OUT_OF_MEMORY` error will be effectively solved by setting this to a higher value, for example 4000MB:

```bash
resources:
  mem_mb=lambda wildcards, attempt: attempt * 4000
```

> [!NOTE]
> Even though you call snakemake using the `--resources mem_mb=100000` flag, the memory might still be limited in the individual rules.
>
> Thus, only changing this flag will likely not fix the issue.

> [!TIP]
> The same logic applies to other resource limitations as for the above mentioned `mem_mb`. When snakemake rules display errors for other
> resource limitations, like `runtime`, you can likely change this in the **resources** directive of the problematic rule.

### Other issues

If there are issues, the documentation of the original in-house workflow might be helpful resolving the issue.