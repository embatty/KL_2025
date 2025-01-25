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




