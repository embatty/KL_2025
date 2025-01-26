# Computation Practical 8: Analysing phylogenetic trees

Module Developers: Dr. Pakorn Aiewsakun and Dr Francesc Coll I Cerezo 

## Table of Contents
1. [Introduction to phylogenetic analysis](#intro)
2. [Expected learning outcomes](#outcomes)
3. [Creating a multiple sequence alignment](#alignment)
4. [Estimating a phylogeny](#estimating)
5. [Phylogeny estimation with bootstrapping](#bootstrapping)
6. [Rooting your tree](#rooting)
7. [How to read a phylogenetic tree?](#readtree)
8. [Relatedness](#relatedness)
9. [Trait evolution](#trait)
10. [Types of phylogenetic groups](#types)
11. [Detecting potential conflicting evolutionary signals within the MSA](#types)
12. [Phylogenetic network reconstruction](#network)
13. [Take-home messages](#messages)
14. [Answers to exercises on interpreting phylogenetic trees](#answers)
15. [Bibliography](#biblio)

---

## Introduction to phylogenetic analysis <a name="intro"></a>

**Phylogenetic analysis** is central to many areas of modern microbiological research, from pathogen classification to examination of pathogen-host co-evolution and tracking of transmission of infectious diseases. The main result from a phylogenetic analysis is a **phylogeny**, also known as **phylogenetic tree**, or simply a tree, depicting evolutionary relationships of a set of taxa – these can be different species or different strains of the same species, as it is often the case in genomic epidemiology studies. **Figure 1** to familiarise with commonly used phylogenetic terminology (e.g., internal nodes, branches, etc,). 

![](images/phy_Figure_1.jpg)  
**Figure 1. Basic phylogenetic terminology. This figure has been reproduced from [Aiewsakun 2024](https://doi.org/10.1016/B978-0-323-99886-4.00013-2).** 

Ideally, we would like to have a complete knowledge of the entire genealogical history of the investigated taxa to draw their phylogeny, but such information is almost always impossible to ascertain. In practice, a phylogeny is commonly estimated from similarity of the investigated organisms’ **orthologous characters** (i.e., biological features with their most recent diversification event coincides with that of the organisms bearing them), with the assumption that the more similar the characters, the more likely that the organisms bearing them would be more closely related (i.e., sharing a more recent common ancestor). Theoretically speaking, any orthologous characters can be used to reconstruct a phylogeny under this framework; however, researchers nowadays most commonly use **DNA sequence data**, such as genes or genome sequences, for this purpose, due to its relative ease of generation, curation, interpretation, transferability, and reproducibility. A phylogeny estimated from molecular sequences is commonly referred to as a **molecular phylogeny**.

In the context of infectious disease epidemiology, a phylogenetic tree is commonly used to depict estimated evolutionary relationships between strains of the same bacterial species, which in turn can be used to track their transmission. Changes of population size can also affect how evolutionary changes accumulate on molecular sequences, and thus a phylogeny can also be analysed to estimate epidemiology dynamics as well. While bacteria typically reproduce its genome with high fidelity (estimated at 1 in 10 million to 1 in a billion base substitutions per nucleotide per generation, very much comparable to that of us human really), random errors (i.e., **mutations**) in DNA replication may still occur, and they can be passed down from one generation to the next. Some mutations are beneficial, while some can be neutral or deleterious. As time goes by, the frequency of beneficial mutations that increase the bacterial fitness (i.e., increase the chance of survival and / or reproductive success of the bacteria in the environment that they find themselves in) will tend to increase in the population, while the frequency of deleterious mutations, i.e., those that decrease the fitness, will tend to decrease as bacteria carrying them will tend to be outcompeted, out-survived, and out-reproduced by others. This in turn will result in continuous change of genetic composition of the bacterial population over time. This process is now known as **evolution by natural selection**, first proposed by Charles Darwin in 1859. Although the mutation rate of bacteria is not high, given their typically short generation times and large population size, a bacterial population will nearly always have mutants around, allowing natural selection to operate rather efficiently. Combined with the usually strong selection pressure that bacteria face constantly, for example from the host immune response, antimicrobial drugs, and other competing organisms, the genetic composition of a bacterial population can and do change substantially within a relatively short amount of time. It is all of these features combined that make molecular phylogenetic analysis of bacteria populations feasible and meaningful.

A bacterial phylogeny is typically estimated from a (set of) orthologous gene(s), more commonly from **polymorphic sites**, i.e., sites showing multiple forms of molecular variants within the population. The number and pattern of shared evolutionary changes between bacterial strains can be used to reconstruct their genealogical and evolutionary relationships. The advancement of sequencing technology nowadays also makes it possible to ‘read’ the entire bacterial genome virtually within a day at an affordable price. This allows researchers to estimate a bacterial phylogeny from genes or polymorphic sites sampled across their whole genome, providing the ultimate level of resolution possible to discriminate between closely related strains, and in turn the finest disease transmission history.

## Expected learning outcomes <a name="outcomes"></a>

After this practical session, you should be able to:
- familiarise with phylogenetic concepts and nomenclature;
- reconstruct a phylogeny from whole-genome sequences using a maximum likelihood framework;
- interpret a phylogeny and assess phylogenetic uncertainty;
- identify and mask recombination from whole-genome sequence alignments;
- understand how to handle conflicting evolutionary signals with phylogenetic network analysis;

## Creating a multiple sequence alignment <a name="alignment"></a>

The very first and arguably the most critical step in reconstructing a phylogeny from molecular sequences, i.e., **molecular phylogenetic reconstruction**, is to align molecular sequences (most often DNA sequences) of the studied organisms together to create a **multiple sequence alignment (MSA)**. With an MSA, we can then estimate the degree of organisms’ (dis)similarity or infer the process of evolutionary changes to estimate their phylogenetic tree. It is very important to note that all positions within the MSA should be **homologous positions** to ensure that we compare ‘like with like’, otherwise the resultant tree will be meaningless. While this might sound simple, in reality, it can oftentimes be difficult, if not entirely impractical, to ensure that all positions within our MSA are homologous positions, especially when the number of analysed organisms is large and they are highly diverse. When this happens, this issue must be dealt with utmost care, making sure that there are as many correctly aligned positions (or conversely as few potentially misaligned positions) as possible in the MSA, in order to obtain the most reliable and meaningful tree. Finding homologous sequences may be especially difficult when dealing with large and highly diverse organisms. In the context of genomic epidemiology investigations **whole-genome sequencing** of multiple strains of the same pathogen often leads to large numbers of **homologous core genes**, i.e., gene sequences that are shared among most strains of the same species and can be used to generate MSA.

**Figure 2** illustrates the common workflow to generate an MSA from a collection of bacterial strains. Generally, bacterial DNA is extracted from a single colony picked from culture plates (therefore commonly referred as to ‘**isolate**’), followed by sequencing library preparation and whole-genome sequencing using rapid benchtop sequencers. Raw sequence data generated by sequencers are then processed using bioinformatic and genomic pipelines, which generally involve read cleaning, and mapping them to a reference genome to reconstruct the isolate’s DNA sequence along the whole bacterial chromosome. With good quality control, mapping the reads of multiple sequenced isolates to the same reference genome is often the way to create an MSA while ensuring that all positions are **homologous positions**. 

![](images/phy_Figure_2.png)  
**Figure 2. Common workflow (bottom) to generate a multiple sequence alignment (top) from a collection of bacterial strains.**

**Polymorphic sites** in the MSA provide information useful for evolutionary relationship inference (i.e., the more the similar the sequences among ‘orthologous’ polymorphic sites, the more likely they share a recent common ancestor), whereas monomorphic sites (nucleotide positions with all individuals showing the same molecular variants) are generally ignored (although it is a good practice to take these into account when specifying equilibrium base frequencies in the phylogenetic estimation, see below). **Figure 3** shows how polymorphic sites may inform phylogenetic estimation.

![](images/phy_Figure_3.png)  
**Figure 3. How polymorphic sites may inform phylogenetic inference.** An example of a simple MSA of eight sites from four strains, which include monomorphic (squared) and polymorphic sites. Genetic changes in the phylogenetic tree are showed as coloured vertical rectangles on the branch where they originated. The identification of genetic changes (alleles) that are unique and common to multiple taxa (strains) are used to group them into phylogenetic clusters in a hierarchical manner with the goal of constructing the most plausible genealogical relationships between strains.

In the previous session, you we run **Snippy** to map the short Illumina reads and **call SNPs** for each individual CPE *Klebsiella pneumoniae* isolate from the suspected hospital outbreak. Now will attempt to construct a molecular phylogeny from this so-called *SNP alignment*.

First, concatenate all Snippy consensus sequences to create a *multi-sequence alignment (MSA)*. Remember all isolates were mapped to the same reference genome (the chromosome of ST78 cpe058 isolate). Here, we will create a *whole-genome alignment* by replacing SNPs and missing alleles called by Snippy along the DNA sequence of the reference genome for each isolate, so all resulting consensus sequences will have the same length (i.e., the length of the chromosome).



