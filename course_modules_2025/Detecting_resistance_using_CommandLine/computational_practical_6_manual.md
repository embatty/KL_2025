# Computational Practical 6: Detecting antimicrobial resistance from bacterial genomes


## Table of Contents
1. [Introduction](#intro)
2. [Expected learning outcomes](#outcomes)
3. [Bacterial strains to be analysed](#strains)
4. [WGS-based prediction of AMR using AMRFinderPlus](#amrfinder)
5. [WGS-based prediction of AMR using ResFinder](#resfinder)
6. [AMR detection using Pathogenwatch](#pathogenwatch)
7. [AMR detection using ResFinder Website – Optional](#resfinderw)
8. [AMR detection using CARD RGI Website – Optional](#card)
9. [Adding genotypic antibiogram of CPE strains into EpiCollect](#epicollect)

---

## Introduction <a name="intro"></a>

Growing rates of antimicrobial resistance make antibiotic susceptibility testing (AST) increasingly needed to ensure the right antibiotics are prescribed for patients with bacterial infections. Determining antibiotic susceptibility is preferred over empiric therapy, wherein typically broad-spectrum drugs are used without a definitive confirmation of the infectious agent and which antibiotics infectious bacteria are resistant to. Data collected on antibiograms (strains’ full susceptibility pattern) can also be used for surveillance purposes and, in turn, inform empiric therapy.

AST is routinely performed using culture-based techniques in clinical diagnostic laboratories, frequently disk diffusion, broth microdilution and gradient diffusion (i.e., E-test). As antibiotic resistance is genetically encoded, i.e. mediated by acquisition of new genes, gene copy number, or mutations in regulatory and coding regions of existing chromosomal genes, molecular tests have been developed to target the detection of such genetic markers. In the last decade, whole-genome sequencing has emerged as an alternative technology to both culture and targeted molecular tests for the detection of AMR as it can, in principle, detect all AMR genetic determinants and predict resistance to all antibiotics in a single experiment. The accuracy of genotypic predictions depends on the availability of: (1) accurate databases of AMR genetic determinants, (2) large collections of whole-genome sequenced strains with AST measurements to assess the diagnostic accuracy of such catalogues, and (3) automated genome analysis and interpretation tools.

Mutational (chromosomal) resistance is the main driver of acquired resistance in certain bacterial species, such as *Mycobacterium tuberculosis* and *Helicobacter pylori*, or for particular antibiotics, especially to synthetic agents such as fluoroquinolones and oxazolidinones. Resistance mutations are vertically transmitted, i.e., via clonal reproduction of bacteria, or can be transmitted horizontally via homologous recombination between different strains. Gene-mediated resistance is the main driver of acquired resistance in certain bacterial species, particularly in gram-negatives. Resistance genes can be horizontally transmitted (via mobile genetic elements such as plasmids) and vertically transmitted via clonal reproduction of bacteria, particularly stable if integrated into the chromosome. In some bacterial species, chromosomal and gene-mediated resistance are equally common (e.g., *Staphylococcus aureus*). Resistance to the same antibiotic can be conferred by both mutations and acquired genes (e.g., fusidic acid in Staphylococcus aureus, colistin resistance in *Escherichia coli*).

Over the years, several global studies have identified the genes and mutations that confer resistance to particular antibiotics. There are several databases such as the [Comprehensive Antimicrobial Resistance Database (CARD)](https://card.mcmaster.ca/), [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/), [AMRFinde](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) or [Pathogenwatch](https://pathogen.watch/) that contain information about the genes and mutations that confer resistance. The use of these databases and tools depends on the species and mechanisms of resistance one is interested in.

## Expected learning outcomes <a name="outcomes"></a>

After this practical session, you should be able to:
- Run broadly used command-line bioinformatic tools like AMRFinder and ResFinder to detect AMR from bacterial genome data;
- Interpret the AMR reports and AMR genetic determinants identified by these tools and their relationship with phenotypic AMR;
- Identify the strengths, limitations, and differences in functionaly and predictions among these tools;
- Employ popular online-based tools like Pathogenwatch or ResFinder to detect AMR from bacterial genome data.

## Bacterial strains to be analysed in this practical <a name="strains"></a>

Table 1 contains the list of strains to be analysed in this practical from the CPE outbreak we will be investigating this year. Table 2 contains additional strains to be analysed (optionally, if time allows) sourced from key studies on the genomic epidemiology of methicillin-resistant *Staphylococcus aureus* (MRSA) ([Holden *et al.* 2013](https://doi.org/10.1101/gr.147710.112)) and extensively drug-resistant (XDR) *Salmonella typhi* ([Klemm *et al.* 2018](https://doi.org/10.1128/mbio.00105-18)). In this practical we will use three different command-based tools (AMRFinder, ResFinder and CARD RGI) to identify AMR genetic determinants from whole-genome sequences. In later computational practicals we will explore the epidemiology of outbreaks and epidemic clones of these bacteria, and study the strains investigated in this practical in a broader context. We will also identify the genomic context of AMR genes, for example, if they are carried on mobile genetic elements.

Table 1 CPE strains to be analysed in this practical
| Species | Study and origin | Strain Id | Illumina accession | Assembly file name |
| :---    | :---             | :---      | :---               | :---      |     
| *K. pneumoniae* | Roberts *et al.* 2024, CPE strain | cpe004 | ERR4095909 | cpe004_Kpn-ST78-NDM1.fasta |
| *E. coli* | Roberts *et al.* 2024, CPE strain | cpe069 | ERR5386320 | cpe069_Eco-NDM1.fasta |

Table 2 Additional strains to be analysed (optional).
| Species	| Study and origin | Strain Id | Genome accession | Assembly file name | 
| :---    | :---             | :---      | :---               | :---      |     
| *S. aureus*	| Holden *et al.* 2013, Berlin (Germany), 2007, ST22 EMRSA-15 | 07-02477 | ERR017261  | ERR017261.assembly.fa |
| *S. aureus*	| Holden *et al.* 2013, UK, 2005, ST22 EMRSA-15	| HO50960412 | HE681097 (GenBank) | HO50960412.fa |
| *S. typhi*	| Klemm *et al.* 2018 (ACT), Pakistan, 2016, 4.3.1 (H58) XDR | BL0006 | ERR2093245	| ERR2093245.assembly.fa |
| *S. typhi* | Klemm *et al.* 2018, Pakistan (2016) – 4.3.1 (H58) pre-XDR	| Pak60168 | ERR2093329	|ERR2093329.assembly.fa |

If you haven’t done so already, clone the course Github directory into your home directory:

```bash
cd ~
git clone https://github.com/WCSCourses/AMR_2025/
```

Create a new directory for this practical named ‘cp6’ and navigate to this directory:
```bash
mkdir ~/course/cp6/
cd ~/course/cp6/
```

And copy the genome assemblies we will use (i.e., those in Tables 1 and 2) into cp6 directory:
```bash
cp ~/AMR_2025/course_data_2025/cp6/complete_assemblies/cpe004_Kpn-ST78-NDM1.fasta ~/course/cp6/
cp ~/AMR_2025/course_data_2025/cp6/complete_assemblies/cpe069_Eco-NDM1.fasta ~/course/cp6/
```

Copy the genome assemblies of the additional strains in Table 2 (the analysis of these strains is optional):
```bash
cp ~/AMR_2025/course_data_2025/cp6/additional_assemblies/HO50960412.fa ~/course/cp6/
cp ~/AMR_2025/course_data_2025/cp6/additional_assemblies/ERR017261.assembly.fa ~/course/cp6/
cp ~/AMR_2025/course_data_2025/cp6/additional_assemblies/ERR2093245.assembly.fa ~/course/cp6/
cp ~/AMR_2025/course_data_2025/cp6/additional_assemblies/ERR2093329.assembly.fa ~/course/cp6/
```

Also, identify and copy the genome assembly of **your assigned CPE strain** (the one on your EpiCollect sheet).

Finally, launch the course Docker image with cp6 directory mounted to it:

```bash
docker run -p 5900:5900 -it --mount type=bind,source=$HOME/course/cp6/,target=/home/data amr:Dockerfile
```

## 4. WGS-based prediction of AMR using AMRFinderPlus <a name="amrfinder"></a>

### Introduction to AMRFinderPlus

To enable accurate assessment of AMR gene content, as part of a multi-agency collaboration, the National Center for Biotechnology Information (NCBI) in the US developed a comprehensive AMR gene database, the Bacterial Antimicrobial Resistance Reference Gene Database, and AMRFinder, an AMR gene detection tool.9 Recently, NCBI released a new version of AMRFinder, known as AMRFinderPlus that, among several new functionalities, has been expanded to detect point mutations in both protein and nucleotide sequences, and taxon-specific analyses that include, or exclude, certain genes and point mutations for specific taxa. (AMRFinderPlus)[https://github.com/ncbi/amr] is available on as a command-line tool only. In this section we will run AMRFinderPlus on the same strain genomes analysed with ResFinder and CARD RGI in previous sections.

### AMRFinderPlus commands

The only required arguments to run AMRFinderPlus are either ```-p <protein_fasta>``` for proteins or ```-n <nucleotide_fasta>``` for nucleotides. Use ```--help``` to see the complete set of options and flags.

```bash
amrfinder --help
```

Use ‘amrfinder -u’ to download and prepare database for AMRFinderPlus:

```bash
amrfinder -u
```

First, a local database of the latest the latest AMR database must be download.

```bash
mkdir amrfinder_db
amrfinder_update -d ./amrfinder_db
```

After making sure the latest AMR database is downloaded, you can run amrfinder on genome assemblies, as showed in the command line below:

```bash
amrfinder -n cpe004_Kpn-ST78-NDM1.fasta -O Klebsiella_pneumoniae -o cpe004_Kpn-ST78-NDM_amrfinder.txt
```

The command above will run amrfinder on the *Klebsiella pneumoniae* strain cpe004 we created an assembly for in previous practicals.

It should take a couple of minutes for this command to finish.

From the command above, note the following chosen options:
- AMRFinder only supports the processing of input nucleotide sequences in FASTA format (with the ```-n/--nucleotide``` option), and not the analysis of raw reads in fastq format. This means that raw reads must be de novo assembled first.
- The option ```-o/--output``` allows you to choose the name of the output file.
- One of the strengths of AMRFinfer is the option ```-O/--organism``` which can be used to get organism-specific results. For those organisms which have been curated, using ```--organism``` will get optimized organism-specific results, and it is therefore recommended. AMRFinderPlus uses the ```--organism``` for screening for point mutations and to filter out genes that are nearly universal in a group and uninformative.

Use ```amrfinder -l``` to list the organism options supported by AMRFinder:

```bash
amrfinder -l
```

You will find taxa like ‘Klebsiella_pneumoniae’, ‘Staphylococcus_aureus’ or ‘Salmonella’ included among the list of supported organisms.

The command below will execute AMRFinder on our CPE *E. coli* strain of interest (Table 1):
```bash
amrfinder -n cpe069_Eco-NDM1.fasta -O Escherichia_coli -o cpe069_Eco-NDM1_amrfinder.txt
```

Now adapt and run the amrfinder command above on your assigned outbreak strain. First, identify and copy the hybrid assembly of your assigned strain into your working directory. Second, make sure to choose the right organism with the parameter ```-O```.

It time allows, come back to this section later to run AMRFinder on the additional strains:
```bash
amrfinder -n HO50960412.fa -O Staphylococcus_aureus -o HO50960412_amrfinder.txt
amrfinder -n ERR017261.assembly.fa -O Staphylococcus_aureus -o ERR017261_amrfinder.txt
amrfinder -n ERR2093245.assembly.fa -O Salmonella -o ERR2093245_amrfinder.txt
amrfinder -n ERR2093329.assembly.fa -O Salmonella -o ERR2093329_amrfinder.txt
```

### Interpreting AMRFinderPlus results

The table below includes a few rows and some of the columns of the AMRFinderPlus output of *Klebsiella pneumoniae* strain cpe004 (file cpe004_Kpn-ST78-NDM_amrfinder.txt).

| Gene symbol | Sequence name | Element subtype | Subclass |
| :---        | :---          | :---            | :---     |  
| aph(3')-VI	| APH(3')-VI family aminoglycoside O-phosphotransferase	| AMR | AMIKACIN/KANAMYCIN |
| aac(6')-Ib-cr5	| fluoroquinolone-acetylating aminoglycoside 6'-N-acetyltransferase AAC(6')-Ib-cr5 | AMR | AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN |
| aac(6')-Ib-cr5 | fluoroquinolone-acetylating aminoglycoside 6'-N-acetyltransferase AAC(6')-Ib-cr5	| AMR | AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN |
| aac(6')-Ib	| AAC(6')-Ib family aminoglycoside 6'-N-acetyltransferase	| AMR	| AMIKACIN/KANAMYCIN/TOBRAMYCIN |

The column ‘Gene symbol’ indicates the genetic determinant (either acquired gene or point mutation) associated with phenotypic resistance, the latter indicated in the column ‘Subclass’.

Based on AMRFinderPlus output files, fill in the tables in the Word document **Summary of genotypic AMR results - CPE strains.docx** to facilitate comparison of WGS-predicted antibiograms between strains.





