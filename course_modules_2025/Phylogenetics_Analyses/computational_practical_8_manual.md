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
**Figure 1.** Basic phylogenetic terminology. This figure has been reproduced from [Aiewsakun 2024](https://doi.org/10.1016/B978-0-323-99886-4.00013-2).







