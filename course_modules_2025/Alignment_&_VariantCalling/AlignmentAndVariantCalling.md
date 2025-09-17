# Computational Practical: Genome Alignment and SNP Calling in Bacterial Genomes
Original version by Dr Arun Gonzales Decano, modified by Liz Batty September 2025

## Table of Contents
1. [Introduction](#intro)
2. [Prerequisites](#prereq)
3. [Data preparation](#dataprep)
4. [Short Read Alignment and SNP Calling with Snippy](#snippy)
5. [Post-processing and Analysis](#postprocessing)

## Introduction  <a name="intro"></a>

Variant calling involves crucial steps in comparative genomics that researchers to identify genetic variations between bacterial strains. Short-read sequencing data (from Illumina) and long-read sequencing data (from PacBio or Oxford Nanopore) have distinct characteristics, each with specific tools and methods for optimal processing. This tutorial provides step-by-step instructions for performing SNP (Single Nucleotide Polymorphism) and other forms of variant calling in bacterial genomes using short-read sequencing data, with an optional exercise looking at long-read data at the end of the practical.

## Prerequisites <a name="prereq"></a>

Before starting, load up the conda environment for this practical:
`conda activate AlignmentAndVariantCalling`

And check you have the following tools available:
  - fastp
  - snippy
  - samtools
  - bcftools
  - fastqc
  - snpeff

## Data Preparation <a name="dataprep"></a>

We will map Illumina short read data from a strain of _Klebsiella pneumoniae_. The data files you will need (read files, and a reference genome) are already downloaded for you and are in the `AlignmentAndVariantCalling` directory wherever you placed the downloaded data.

<details>
    <summary>Downloading data</summary>
    The following commands were used to download the data directly from the Sequence Read Archive.

    ```
    wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095905/ERR4095905_1.fastq.gz
    wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095905/ERR4095905_2.fastq.gz
    wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095885/ERR4095885_1.fastq.gz
    wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095885/ERR4095885_2.fastq.gz
    ```

    Download a Klebsiella pneumoniae reference genome from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000240185.1/

</details>


## Short Read Alignment and SNP Calling with Snippy <a name="snippy"></a>
First we prepare our data using `fastp` to remove adapters and low-quality bases, as we practiced in the Data QC module.

     ```
     fastp -i ERR4095905_1.fastq.gz -I ERR4095905_2.fastq.gz -o out.ERR4095905_1.fastq.gz -O out.ERR4095905_2.fastq.gz
     ```

[Snippy](https://github.com/tseemann/snippy) is an all-in-one tool for bacterial SNP calling using short-read data. It aligns the reads to a reference genome and calls variants.

Run Snippy and examine the output.
     ```
     snippy --cpus 1 --outdir ERR4095905_snippy --reference cpe058_Kpn-ST78-NDM1.chr.fasta --R1 out.ERR4095905_1.fastq.gz --R2 out.ERR4095905_2.fastq.gz
     ```
This will take about ten minutes to run. All the output files from your Snippy run will be placed in the output directory `ERR4095905_snippy`. Look at the output files and refer to the Snippy homepage to tell you what is in each output file.

 A file ending in VCF is in the [Variant Call Format](https://en.wikipedia.org/wiki/Variant_Call_Format), which is a specific file format designed to hold information on genetic variants. A BAM file and the associated BAI index file are a specific file format which contains the reads aligned to the reference genome. A BAM file is a binary file which must be viewed using the `samtools` application.
     | Filename | Description |
     | snps.vcf | The called SNPs in VCF format |
     | snps.tab | A tabular summary of SNPs |
     | alignment.bam | The aligned reads in BAM format |
     | alignment.bam.bai | The BAM index file |

Now we run snippy-core and examine the output. Snippy-core compares the output of SNP calling on multiple samples run using Snippy to look at the 'core' genome, where we can call variants in all samples. Three previously run on Snippy can be found in the `AlignmentAndVariantCalling/snippy` directory. First, navigate to this directory, and then run Snippy-core across all samples:

  ```
  snippy-core --ref Kpne_HS11286.fna snippy/ERR4095905_snippy/ snippy/ERR4095977_snippy/ snippy/ERR9419473_snippy/
  ```

  Note the formats and contents of each output file. By default all the files produced by snippy-core will starting with `core.`


## Post-processing and Analysis <a name="postprocess"></a>

### 1. Filter SNPs:
We can also apply filters to the VCF files to remove or tag low-quality SNPs using `bcftools`. For instance, to remove SNPs with a depth of coverage (DP) below 30:
     ```
     bcftools filter -e 'FMT/DP<30' snps.vcf > snps.depthfilter.vcf
     ```

How many SNPs were removed through filtering?
<detail><summary>Answer</summary>
You can count the number of lines in the SNP file before and after filtering using the `wc -l` command.
</detail>

### 2. Annotation of SNPs:
Use [snpEff](https://pcingola.github.io/SnpEff/snpeff/running/) to annotate the SNPs.

     ```
     #short reads
      java -Xmx8g -jar snpEff.jar reference filtered_snps_short.vcf > annotated_snps_short.vcf
     ```

### 3. Visualization:
   - Visualize the alignment and SNPs using tools using [Integrative Genomics Viewer (IGV)](https://igv.org/doc/desktop/#DownloadPage/). This cannot be installed using Conda but you can download it directly to your computer.

<details>
    <summary>Explanation of each metric in the pop-up overview of a SNP</summary>

*Basic Information*

ID: This field is empty, indicating no specific identifier (such as a dbSNP ID) is assigned to this SNP.

Chr: NC_009648.1 is the chromosome or contig reference name from the genome assembly.

Position: 2,803,960 is the position on the chromosome where this SNP is located.

Reference: C* indicates the reference allele (the allele found in the reference genome) at this position is "C".

Alternate: T is the alternate allele observed in this SNP, meaning a "C" in the reference is replaced by "T" in the observed variant.

Qual (Quality): 882.977 is a quality score for the variant call. It indicates the confidence in the variant being true; higher scores mean greater confidence.

Type: SNP specifies that this variant is a single nucleotide polymorphism.

Is Filtered Out: No indicates that this SNP has passed any filters applied during variant calling, suggesting it is considered a reliable call.

*Alleles Information*

Alternate Alleles: T is the alternate allele for this SNP.

*Variant Attributes*

QA (Quality of Alternate): 1039 is the sum of base quality scores for the reads supporting the alternate allele (T). Higher values suggest more confidence in the variant call.

AB (Allele Balance): 0 represents the proportion of reads supporting the alternate allele relative to the total depth. Since it is 0, this may mean that the calculation could not be performed, or the data isn't available.

QR (Quality of Reference): 32 is the sum of base quality scores for the reads supporting the reference allele (C). Lower values suggest fewer or less confident reads supporting the reference allele.

Depth: 30 is the total read depth at this position, representing the number of reads covering the SNP site.

RO (Reference Observations): 1 indicates the number of reads supporting the reference allele (C).

TYPE: snp reaffirms that this variant is a single nucleotide polymorphism.

AO (Alternate Observations): 29 is the number of reads supporting the alternate allele (T).
</summary>
</details>

### Long read variant calling
Variant calling from long reads can be performed in a similar way to short reads but using different tools designed for long read sequencing.
<details><summary>Long read variant calling protocol</summary>

Long reads go here
</details>
## References

1. Bush SJ, Foster D, Eyre DW, Clark EL, De Maio N, Shaw LP, Stoesser N, Peto TEA, Crook DW, Walker AS, Wilson DJ. (2020). Genomic diversity affects the accuracy of bacterial SNP calling pipelines. Genome Biology, 21, 20. https://doi.org/10.1186/s13059-019-1921-9

2. Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191

3. Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, Newburger D, Dijamco J, Nguyen N, Afshar PT, Gross SS, Dorfman L, McLean CY, DePristo MA. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987. https://doi.org/10.1038/nbt.4235

4. Cingolani P, Platts A, Wang LL, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80-92. https://doi.org/10.4161/fly.19695

5. Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G, Mesirov JP. (2011). Integrative genomics viewer. Nature Biotechnology, 29(1), 24-26. https://doi.org/10.1038/nbt.1754
