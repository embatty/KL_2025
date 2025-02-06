# Tutorial: Metagenomic Analysis of Bacterial Genomes Using Illumina Short-Read Data (Antimicrobial Resistance in Bacterial Pathogens - Asia 2025)

**Aarthi Ravikrishnan, Arun Decano**

### Learning Objectives

This tutorial provides a step-by-step command-line workflow for analyzing bacterial genomes from Illumina short-read metagenomic data. By following these steps, you’ll be able to:

1. Conduct quality control on Illumina reads.
2. Assemble short reads into contigs using **metaSPAdes**.
3. Bin contigs into genomes (Metagenome-Assembled Genomes or MAGs).
4. Taxonomically classify and annotate these genomes.
5. Assess antimicrobial resistance (AMR) potential and visualize taxonomic data.
6. Understand how to perform high-throughput analyses using a pre-made Snakemake workflow (https://github.com/aarthi31/amr-course-new).

### Prerequisites

To complete this tutorial, you need:

- **Basic Linux command-line** knowledge.
- **Installed tools**: FastQC, Fastp, seqtk, Bbmap, metaSPAdes, MetaBAT2, Kraken2, CheckM, Prokka, Abricate (for AMR prediction), and Pavian or Krona for visualization.

### Dataset

Raw files were grabbed from Guo, X., Tang, N., Lei, H., Fang, Q., Liu, L., Zhou, Q., & Song, C. (2021). Metagenomic Analysis of Antibiotic Resistance Genes in Untreated Wastewater From Three Different Hospitals. Frontiers in microbiology, 12, 709051. [https://doi.org/10.3389/fmicb.2021.709051](https://doi.org/10.3389/fmicb.2021.709051).

**Command to obtain raw sequences**

```bash
# Raw
wget <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR142/072/SRR14297772/SRR14297772_1_ds.fastq.gz>
wget <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR142/072/SRR14297772/SRR14297772_2.fastq.gz>

```

For this workshop, we will be using a custom set of sequences (which were spiked-in) that can be obtained from here ([https://tinyurl.com/mvr3d263](https://tinyurl.com/mvr3d263)).

You can also find these sequences within the `metagenomics/mgx_sequences` folder.

```
# Spiked: c/o Aarthi and Arun: Download from <https://tinyurl.com/mvr3d263>
# Filenames: SRR14297772_cpe107_1_ds.fastq.gz and SRR14297772_cpe107_2.fastq.gz

```
Copy these sequences to the current working directory `metagenomics`.

```
cp mgx_sequences/SRR14297772_cpe107_* .
```

## Step 1: Quality Control of Illumina Reads

Quality control ensures your Illumina reads are suitable for assembly. **FastQC** will identify quality issues, and **Fastp** will trim low-quality bases and adapters.

1. **Run FastQC to Assess Read Quality**:
   First make an output directory wihtin `metagenomics` folder.
   ```
    mkdir fastqc_output/
   ```
   Next run fastqc on the sequences in this directory. Note the use of '*' as a wildcard.  
    ```
    fastqc *.fastq.gz -o fastqc_output
    
    ```
    
    - **What It Does**: FastQC generates quality metrics and HTML reports, including GC content, read length distribution, and base quality scores.
    - **Output**: HTML reports in the `fastqc_output` directory. Open these to check for any quality issues.
3. **Clean Reads Using Fastp and remove host reads using Hostile**:

This is a general command used to run fastp. Replace the input and output files based on the files you have.

**General command**

```bash
fastp -i raw_reads/sample_R1.fastq -I raw_reads/sample_R2.fastq -o trimmed_sample_R1.fastq -O trimmed_sample_R2.fastq -h fastp_report.html -j fastp_report.json --length_required 50
```

For the dataset given to you, try the following command to run `fastp`.

```bash
# We will use the spiked pair
fastp -i SRR14297772_cpe107_1.fastq.gz -I SRR14297772_cpe107_2.fastq.gz -o SRR14297772_cpe107_1_ds_filtered.fastq.gz -O SRR14297772_cpe107_2_ds_filtered.fastq.gz -h fastp_report.html -j fastp_report.json --length_required 50

```

- **What It Does**: Fastp removes adapters, trims low-quality bases, and discards reads shorter than 50 bp. Hostile gets rid of host (mainly) human reads.
- **Key Options for Fastp**:
    - `h` and `j`: Generate HTML and JSON reports with trimming metrics.
    - `-length_required`: Discards reads shorter than 50 bp.
- **Hostile**: [https://github.com/bede/hostile](https://github.com/bede/hostile)

```bash
# Create and activate a conda env 
conda create -y -n hostile -c conda-forge -c bioconda hostile 
conda activate hostile
conda activate --stack amr # to access packages from the amr env 
 # Run Hostile on paired short reads 
hostile clean --fastq1 SRR14297772_cpe107_1_ds_filtered.fastq.gz --fastq2 SRR14297772_cpe107_2_ds_filtered.fastq.gz -o - > SRR14297772_cpe107.interleaved.fastq
 
```

```bash
# Bin interleaved fastq files into clean.fastq1 and clean.fastq2 using seqtk
seqtk seq -1 SRR14297772_cpe107.interleaved.fastq > clean.SRR14297772_cpe107_1.fastq
seqtk seq -2 SRR14297772_cpe107.interleaved.fastq > clean.SRR14297772_cpe107_2.fastq

#Compress all fastq files (pigz offers fast compression. You can also use gzip)
pigz SRR14297772_cpe107.interleaved.fastq
pigz clean.SRR14297772_cpe107_1.fastq
pigz clean.SRR14297772_cpe107_2.fastq

```

**IMPORTANT**
In this tutorial, we will be using the files from cleaned from fastp. However, in the practical applications, please ensure that you have removed the host reads to ensure your downstream analyses is on microbial reads only.

1. **Verify Trimming Results**:
Rerun FastQC and run Multiqc on clean fastq files to confirm improvements.

## Step 2: Downsample Reads (Optional)

Sometimes this step is done -- but it is not mandatory.
For this tutorial, we will skip this step.

Reduce the number of reads in the dataset while preserving the diversity of the sample.
Use **seqtk** sample for random subsampling of reads to a desired percentage or absolute number.

```bash
seqtk sample -s100 input.fastq 0.1 > downsampled.fastq

#-s100 specifies a random seed for reproducibility.
#0.1 specifies 10% of the total reads (adjust according to needs).

```

Use [**bbnorm.sh**](http://bbnorm.sh/) in BBMap for read normalization, which reduces redundancy while keeping unique reads.

```bash
bbnorm.sh in=input.fastq out=downsampled.fastq target=20 min=2

# target=20 controls coverage normalization depth (can be adjusted based on data).

```

## Step 3: Metagenome Assembly with metaSPAdes

**metaSPAdes** is optimized for metagenomic data and assembles reads into contigs, reconstructing genome fragments from complex microbial communities.

1. **Run metaSPAdes**:

**If you have removed host reads**

```bash
metaspades.py -1 clean.SRR14297772_cpe107_1.fastq.gz -2 clean.SRR14297772_cpe107_2.fastq.gz -o clean_SRR14297772_cpe107_metaspades_output/ --only-assembler # fastp and hostile-cleaned

```

**With output from fastp**

```bash
metaspades.py -1 SRR14297772_cpe107_1_ds_filtered.fastq.gz -2 SRR14297772_cpe107_2_ds_filtered.fastq.gz -o SRR14297772_cpe107_metaspades_output/ --only-assembler # fastp-cleaned only

```

- **What It Does**: metaSPAdes assembles contigs by building a de Bruijn graph adapted for metagenomic data.
- **Output**: Assembled contigs are saved in the `metaspades_output/` directory.
- **Warning**: This step is usually time consuming and memory intensive.

**A**. **Review Assembly Results**:
Inspect `contigs.fasta` in `metaspades_output/` to check contig lengths and quality. What are the N50 and L50 values? [Hint: use python script `find_assembly_stats.py` from the `Metagenomics` folder]. 

_Other tools_ : Check [Quast](https://github.com/ablab/quast).

## Step 4: Binning the Contigs

Binning groups contigs into bins representing putative genomes. **MetaBAT2** performs binning based on sequence composition and read coverage.

1. **Run MetaBAT2**:
    
    ```
    metabat2 -i metaspades_output/contigs.fasta -o bins_folder/bin -m 1500
    
    ```
    
    - **What It Does**: MetaBAT2 clusters contigs into bins that represent draft genomes.
    - **Key Option**:
        - `m 1500`: Sets the minimum contig length to 1500 bp for binning.
2. **Examine Binning Results**:
The binned genomes are saved in `bins_folder/`, with each bin corresponding to a draft genome.

## Step 5: Quality Assessment of Bins

Use **CheckM** to evaluate the quality of binned genomes, assessing completeness and contamination based on conserved marker genes.

1. **Run CheckM**:
    
    ```
    checkm lineage_wf -x fa -t 8 bins_folder/ checkm_output/ --pplacer_threads 8
    
    ```
    
    - **What It Does**: CheckM evaluates each bin for genome completeness and contamination.
    - **Interpret Results**: Bins with >90% completeness and <5% contamination are considered high-quality.
  
Q: What is the completeness and contamination of bins? [**Hint**: Look at bin_stats_ext.tsv]

## Step 6: Taxonomic Classification of Contigs (MAGs) and/or Reads

Classify contigs to identify their taxonomic origin with [**Kraken2**](https://ccb.jhu.edu/software/kraken2/), which compares contigs to a taxonomic database.

1. **Run Kraken2 to identify microbial diversity**:
    
    ```bash
    kraken2 --db /home/data/kraken2/ --threads 8 --output kraken_output.txt --report kraken_report.txt metaspades/contigs.fasta
    ```
    
    - **What It Does**: Kraken2 assigns taxonomic classifications by matching sequences against a reference database.
    - **Output**: Results are saved in `kraken_output.txt` with a summary report in `kraken_report.txt`.
    - Alternatively, run Kraken 2 on the clean reads.
        
        ```
        kraken2 --db /home/data/kraken2/ clean.SRR14297772_cpe107_1.fastq.gz clean.SRR14297772_cpe107_2.fastq.gz --threads 8 --output kraken_output_reads.txt --report kraken_report_reads.txt
        
        ```
        

## Step 7: Genome Annotation

Annotate each bin to identify genes and other genomic features using **Prokka**.

1. **Run Prokka for Genome Annotation**:
    
    ```
    prokka --outdir annotation_output --prefix bin_1_annotation bins_folder/bin.1.fa
    
    ```
    
    - **What It Does**: Prokka annotates genes and functional elements in each bin.
    - **Output**: Annotations are saved in `annotation_output/`.

## Step 8: AMR Prediction

Identify antimicrobial resistance (AMR) genes using **ABRicate**, which screens genomes against known AMR gene databases. By default it uses NCBI database, which is a subset of the AMRFinderPlus database to do AMR gene detection. To exploit the complete functionality of AMR prediction, use AMRFinderPlus. See note [here](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/). 

1. **Run ABRicate for AMR Prediction**:

   List all available databases:
   ```
   abricate --list
   ```

   To run ABRicate on the bins [Default database is ncbi]
   ```
   abricate bins_folder/bin.1.fa > abricate_output.txt
   ```
   - **What It Does**: ABRicate searches for known AMR genes by comparing genome sequences against the ResFinder database.
   - **Output**: Results in `abricate_output.txt` list detected AMR genes, their identities, and resistance classes.
   
**Question**: Can you compare the results from Prokka and ABRicate? Do you find anything in common? Can you explain your observation?

## Step 9: Visualization of Taxonomy with Pavian or Krona

Visualize taxonomic classifications interactively using **Pavian** or **Krona**:

- **Pavian** provides an interactive web-based interface (based on R).
- **Krona** produces circular, hierarchical plots for exploring multi-level taxonomic data.

### Option 1: Visualization with Pavian

1. **Set Up Pavian for Visualization**:
    
    ```
    pavian server
    # or upload the kraken reports on the web server: [<https://shiny.hiplot.cn/pavian/>](https://fbreitwieser.shinyapps.io/pavian/)
    
    ```
    
    - **What It Does**: Pavian launches a local server for visualizing taxonomic classifications in your web browser.
    - **Upload**: Load `kraken_report.txt` into Pavian for an interactive visualization.

### Option 2: Visualization with Krona
Documentation: https://github.com/marbl/Krona/wiki/Installing
```
# Installation
git clone https://github.com/marbl/Krona.git
cd Krona/KronaTools
./install.pl
```
1. **Convert Kraken2 Output for Krona**:
    
    ```
    cut -f2,3 kraken_output.txt > krona_input.txt
    
    ```
    
    - **What It Does**: Extracts only the taxonomic ID and classification from Kraken2 output, creating a format compatible with Krona.
2. **Generate Krona Plot**:
    
    ```
    # Update taxonomy index
    ktUpdateTaxonomy.sh
    
    ```
    
    ```
    ktImportTaxonomy -m 1 -o krona-test.html
    ktImportTaxonomy krona_input.txt -o krona_output.html
    
    ```
    
    - **What It Does**: Krona generates an HTML file (`krona_output.html`) with a multi-level circular plot for exploring taxonomic data.
    - **View the Plot**: Open `krona_output.html` in any web browser.
3. **Exploring the Krona Plot**:
    - Click on different sections to zoom into taxonomic levels.
    - Hover over sections to view specific taxonomic details.

### Exercises

1. Test this workflow on the sample Illumina dataset we provided.
2. Experiment with different assembly parameters in metaSPAdes.
3. Compare taxonomic classifications generated by Kraken2 with other tools, such as Centrifuge.
4. Visualize different datasets in Pavian or Krona to identify any microbial community trends or outliers.
