# Genome Alignment and SNP Calling in Bacterial Genomes: A Tutorial for AMR Analysis Workshop in Bangkok, Thailand (February 2025)
Dr. Arun Gonzales Decano

## Introduction

Variant calling involves crucial steps in comparative genomics that researchers to identify genetic variations between bacterial strains. Short-read sequencing data (from Illumina) and long-read sequencing data (from PacBio or Oxford Nanopore) have distinct characteristics, each with specific tools and methods for optimal processing. This tutorial provides step-by-step instructions for performing SNP (Single Nucleotide Polymorphism) and other forms of variant calling in bacterial genomes using both short-read and long-read (nanopore) sequencing data.

## Prerequisites

Before starting, ensure you have the following software installed:
- For short reads:
  - fastp
  - snippy
  - samtools
  - bcftools
- For long reads:
  - minimap2
  - samtools
  - Medaka
  - DeepVariant
- General tools:
  - FastQC for quality control
  - MultiQC for aggregating QC reports
  - snpEff

## Data Preparation

Grab the raw fastq files from ENA. We will map Illumina short read data associated with the outbreak case isolate you've been assigned. If you have been assigned a case isolate of a different bacterial species other than Klebsiella pneumoniae, use the following run accessions as test samples to execute the commands below: ERR4095905 and ERR4095885. Use the E. coli sample, ERR3284704 for the long read mapping and variant calling and EC958 as the reference genome.


```
# Short reads:
cd /home/data/VC_module/
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095905/ERR4095905_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095905/ERR4095905_2.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095885/ERR4095885_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR409/005/ERR4095885/ERR4095885_2.fastq.gz

# Long reads
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR328/004/ERR3284704/ERR3284704.fastq.gz
# Download the E. coli reference from https://www.ncbi.nlm.nih.gov/nuccore/NZ_HG941718.1?report=fasta and rename with Ecoli_ref.fasta.
mv sequence.fasta Ecoli_ref.fasta
```

### 1. Quality Control:
   - Use FastQC to assess the quality of your sequencing reads.
     ```
     fastqc *.fastq.gz
     ```
   - Aggregate the results with MultiQC.
     ```
     multiqc .
     ```

### 2. Cleaning (if necessary):
   - Use Fastp for short reads to remove adapters and low-quality bases.
     
     ```
     fastp -i ERR4095905_1.fastq.gz -I ERR4095905_2.fastq.gz -o out.ERR4095905_1.fastq.gz -O out.ERR4095905_2.fastq.gz
     fastp -i ERR4095885_1.fastq.gz -I ERR4095885_2.fastq.gz -o out.ERR4095885_1.fastq.gz -O out.ERR4095885_2.fastq.gz
     ```
### 3. Make sure you're in the correct working directory where the files are accessible. 

 
     cd /home/data/VC_module/ 
     # /home/data/VC_module/Practical_outputs/ contains output files and subfolders from pre-runs.
     

    
## Short Read Alignment and SNP Calling with Snippy

Snippy is an all-in-one tool for bacterial SNP calling using short-read data. It aligns the reads to a reference genome and calls variants.
Documentation: https://github.com/tseemann/snippy

### 1. Run Snippy:
   - Use snippy to perform alignment and SNP calling in one step.
     ```
     #PSA
     snippy --cpus 4 --outdir snippy_output --reference reference.fasta --R1 output_R1_paired.fastq --R2 output_R2_paired.fastq
     snippy --cpus 4 --outdir ERR4095905_snippy --reference cpe058_Kpn-ST78-NDM1.chr.fasta --R1 out.ERR4095905_1.fastq.gz --R2 out.ERR4095905_2.fastq.gz
     snippy --cpus 4 --outdir ERR4095885_snippy --reference cpe058_Kpn-ST78-NDM1.chr.fasta --R1 out.ERR4095885_1.fastq.gz --R2 out.ERR4095885_2.fastq.gz
     ```
     ```
     #MSA
     snippy-core --prefix core --ref cpe058_Kpn-ST78-NDM1.chr.fasta *_snippy # or specify the snippy folders: ERR4095885_snippy ERR4095905_snippy
     ```

### 2. Examine the Output: Note the formats and contents of each output file. Compare results, identify key features, and observe any patterns or unexpected findings.
     
   - The output directory will contain various files, including:
     - snps.vcf: The called SNPs in VCF format.
     - snps.tab: A tabular summary of SNPs.
     - alignment.bam: The aligned reads in BAM format.
     - alignment.bam.bai: The BAM index file.
     
     ```
     # PSA
     cd ERR4095905_snippy/
     ls
     ```
     ```
     #MSA
     ls core.*
     #core.aln core.tab core.tab core.txt core.vcf
     grep ">" core.full.aln # view names of aligned genomes
     ```

## Long Read Alignment with minimap2 and SNP Calling with Medaka or DeepVariant.

### Run Medaka:
  - Use Medaka to call SNPs from the long-read BAM file.
    Documentation: https://github.com/nanoporetech/medaka
  
    ```
    #medaka_variant -i longread.input.fastq.gz -r reference.fasta
    medaka_variant -i ERR3284704.fastq.gz -r Ecoli_ref.fasta
    
    ```

### Run DeepVariant: (Optional)
   - Use minimap2 for aligning long reads.
     ```
     minimap2 -a reference.fasta long_reads.fastq > aligned_long_reads.sam
     ```
     
   - Convert and sort the alignment.
     ```
     samtools view -S -b aligned_long_reads.sam > aligned_long_reads.bam
     samtools sort aligned_long_reads.bam -o sorted_long_reads.bam
     ```
     
   - Index the BAM file.
     ```
     samtools index sorted_long_reads.bam
     ```
   - Remove duplicates.
     ```
     samtools rmdup sorted_long_reads.bam deduplicated.bam
     ```

   - Use DeepVariant to call SNPs from the long-read BAM file.
     ```
     # Assuming DeepVariant is installed and properly set up

     run_deepvariant --model_type PACBIO --ref reference.fasta --reads sorted_long_reads.bam --output_vcf dv_output.vcf --output_gvcf dv_output.g.vcf --num_shards 4
     ```
### Use PEPPER DeepVariant for nanopore reads: https://github.com/kishwarshafin/pepper  (Optional)

```
## Pull the docker image.
sudo docker pull kishwars/pepper_deepvariant:r0.8

# Run PEPPER-Margin-DeepVariant
sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont_r10_q20
```


## Post-processing and Analysis  (Optional)

### 1. Filter SNPs:
   - Apply filters to the VCF files to ensure high-quality SNPs using bcftools. Remove SNPs/Indels with MQ<30 and DP<10.
     ```
     bcftools filter -s LowQual -e '%QUAL<30 || DP<10' snippy_output/snps.vcf > filtered_snps_short.vcf
     bcftools filter -s LowQual -e '%QUAL<30 || DP<10' medaka.annotated.vcf > filtered_snps_long.vcf
     ```

### 2. Annotation of SNPs:
   - Use snpEff to annotate the SNPs: https://pcingola.github.io/SnpEff/snpeff/running/
     
     ```
     #short reads
      java -Xmx8g -jar snpEff.jar <reference_database> filtered_snps_short.vcf > annotated_snps_short.vcf

     #long reads
      java -Xmx8g -jar snpEff.jar reference filtered_snps_long.vcf > annotated_snps_long.vcf
     ```

### 3. Visualization:
   - Visualize the alignment and SNPs using tools using Integrative Genomics Viewer (IGV).
   - Download IGV from https://igv.org/doc/desktop/#DownloadPage/.
   - Example pop-up info of a selected SNP:
     
ID: .

Chr: NC_009648.1

Position: 2,803,960

Reference: C*

Alternate: T

Qual: 882.977

Type: SNP

Is Filtered Out: No

Alleles:

Alternate Alleles: T

*Variant Attributes*

QA: 1039

AB: 0

QR: 32

Depth: 30

RO: 1

TYPE: snp

AO: 29

   - Explanation of each metric:

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


## References

1. Bush SJ, Foster D, Eyre DW, Clark EL, De Maio N, Shaw LP, Stoesser N, Peto TEA, Crook DW, Walker AS, Wilson DJ. (2020). Genomic diversity affects the accuracy of bacterial SNP calling pipelines. Genome Biology, 21, 20. https://doi.org/10.1186/s13059-019-1921-9

2. Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191

3. Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, Newburger D, Dijamco J, Nguyen N, Afshar PT, Gross SS, Dorfman L, McLean CY, DePristo MA. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987. https://doi.org/10.1038/nbt.4235

4. Cingolani P, Platts A, Wang LL, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly, 6(2), 80-92. https://doi.org/10.4161/fly.19695

5. Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G, Mesirov JP. (2011). Integrative genomics viewer. Nature Biotechnology, 29(1), 24-26. https://doi.org/10.1038/nbt.1754

6. Andrews S. (2010). FastQC: A quality control tool for high throughput sequence data. Available online: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

7. Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354
