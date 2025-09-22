# Installation of ResFinder
You can install ResFinder using the `conda install resfinder` command. It has certain requirements which are difficult
Execute the command below to display all ResFinder arguments and options:
```bash
python -m resfinder --help
```

Next, if not already available, download a local copy of the latest ResFinder databases:
```bash
git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/

# Local databases need to be indexed using kma:
cd resfinder_db
python3 INSTALL.py kma_index
cd ..
cd pointfinder_db
python3 INSTALL.py kma_index
cd ..
cd disinfinder_db
python3 INSTALL.py kma_index
cd ..
```

Set approximate environment bash variables for ResFinder executable to locate these databases.
```bash
export CGE_RESFINDER_RESGENE_DB="/home/data/cp6/resfinder_db";
export CGE_RESFINDER_RESPOINT_DB="/home/data/cp6/pointfinder_db";
export CGE_DISINFINDER_DB="/home/data/cp6/disinfinder_db";
```

Remember to set these variables in any new terminal window. Otherwise ResFinder will exist with the error: ‘Could not locate ResFinder database path’.

Now everything is set to run ResFinder on your terminal screen as shown in the commands below:
```bash
python -m resfinder -ifa cpe004_Kpn-ST78-NDM1.fasta -s "Klebsiella" --acquired --point --outputPath cpe004_Kpn_resfinder
python -m resfinder -ifa cpe069_Eco-NDM1.fasta -s "Escherichia coli" --acquired --point --outputPath cpe069_Eco_resfinder
```
 <!---

Execute the command below to display all ResFinder arguments and options:
```bash
python -m resfinder --help
```

Next, if not already available, download a local copy of the latest ResFinder databases:
```bash
git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/

# Local databases need to be indexed using kma:
cd resfinder_db
python3 INSTALL.py kma_index
cd ..
cd pointfinder_db
python3 INSTALL.py kma_index
cd ..
cd disinfinder_db
python3 INSTALL.py kma_index
cd ..
```

Set approximate environment bash variables for ResFinder executable to locate these databases.
```bash
export CGE_RESFINDER_RESGENE_DB="/home/data/cp6/resfinder_db";
export CGE_RESFINDER_RESPOINT_DB="/home/data/cp6/pointfinder_db";
export CGE_DISINFINDER_DB="/home/data/cp6/disinfinder_db";
```

Remember to set these variables in any new terminal window. Otherwise ResFinder will exist with the error: ‘Could not locate ResFinder database path’.

Now everything is set to run ResFinder on your terminal screen as shown in the commands below:
```bash
python -m resfinder -ifa cpe004_Kpn-ST78-NDM1.fasta -s "Klebsiella" --acquired --point --outputPath cpe004_Kpn_resfinder
python -m resfinder -ifa cpe069_Eco-NDM1.fasta -s "Escherichia coli" --acquired --point --outputPath cpe069_Eco_resfinder
```

IMPORTANT NOTE: if ResFinder database could not be found (```‘Could not locate ResFinder database path’```, you can use the option ```--db_path_res``` to indicate where the directory of such database is:

```bash
python -m resfinder -ifa cpe004_Kpn-ST78-NDM1.fasta -s "Klebsiella" --acquired --point --outputPath cpe004_Kpn_resfinder --db_path_res ./resfinder_db
```

```bash
python -m resfinder -ifa cpe069_Eco-NDM1.fasta -s "Escherichia coli" --acquired --point --outputPath cpe069_Eco_resfinder --db_path_res ./resfinder_db
```
--->

The command line above was used to run ResFinder on the genome assembly of *Klebsiella pneumoniae* cpe004 and *Escherichia coli* cpe069 strains (Table 1). Note the following parameters:
- the option ```-ifa``` is used to indicate that the input genome is provided in FASTA format, following by the path to the genome assembly file we want to analyse;
- the option ```-s``` is used to indicate the bacterial species in the same. This is important for ResFinder to use the antimicrobial panel specific to each bacterial species;
- the option ```--acquired``` is chosen to detected acquired resistance genes, and;
- the option ```--point``` to scan for AMR chromosomal mutations;
- the option ```--outputPath``` allows you to specify the name of the output directory where ResFinder files will be stored

Next, we will run ResFinder on the raw sequencing reads of strain cpe004 (Illumina accession number ERR4095909). These files have already been downloaded, but could be obtained directly from the [ENA](https://www.ebi.ac.uk/ena/browser/view/ERR4095909).

In the ResFinder command below note we used the same options as for sample cpe004 except for ‘-ifq’, used here to specify input fastq file(s). ResFinder assumes the input to be single-end fastq if only one file is provided after ‘-ifq’, and to be paired-end data if two files are provided instead.

```bash
python -m resfinder -ifq ERR4095909_1.fastq.gz ERR4095909_2.fastq.gz -s "Klebsiella" --acquired --point --outputPath cpe004_ERR4095909_resfinder
```

### Optional - if time allows

You can use the commands below to run ResFinder for the *S. aureus* samples:
```bash
python -m resfinder -ifa HO50960412.fa -s "Staphylococcus aureus" --acquired --point --outputPath HO50960412_resfinder
python -m resfinder -ifa ERR017261.assembly.fa -s "Staphylococcus aureus" --acquired --point --outputPath ERR017261_resfinder
```

You can use the commands below to run ResFinder for the *S. typhi* samples:
```bash
python -m resfinder -ifa ERR2093245.assembly.fa -s "Salmonella enterica" --acquired --point --outputPath ERR2093245_resfinder
python -m resfinder -ifa ERR2093329.assembly.fa -s "Salmonella enterica" --acquired --point --outputPath ERR2093329_resfinder
```

Let's now change the delimiter of ResFinder output files to make it easier to open in Excel:
```bash
cat ./HO50960412_resfinder/ResFinder_results_tab.txt | tr '\t' ',' > ./HO50960412_resfinder/ResFinder_results_tab.csv
cat ./HO50960412_resfinder/pheno_table_staphylococcus_aureus.txt | tr '\t' ',' > ./HO50960412_resfinder/pheno_table_staphylococcus_aureus.csv

cat ./ERR017261_resfinder/ResFinder_results_tab.txt | tr '\t' ',' > ./ERR017261_resfinder/ResFinder_results_tab.csv
cat ./ERR017261_resfinder/pheno_table_staphylococcus_aureus.txt | tr '\t' ',' > ./ERR017261_resfinder/pheno_table_staphylococcus_aureus.csv

cat ./ERR2093245_resfinder/ResFinder_results_tab.txt | tr '\t' ',' > ./ERR2093245_resfinder/ResFinder_results_tab.csv
cat ./ERR2093245_resfinder/pheno_table_salmonella_enterica.txt | tr '\t' ',' > ./ERR2093245_resfinder/pheno_table_salmonella_enterica.csv

cat ./ERR2093329_resfinder/ResFinder_results_tab.txt | tr '\t' ',' > ./ERR2093329_resfinder/ResFinder_results_tab.csv
cat ./ERR2093329_resfinder/pheno_table_salmonella_enterica.txt | tr '\t' ',' > ./ERR2093329_resfinder/pheno_table_salmonella_enterica.csv
```

### Interpreting ResFinder results

Based on the ResFinder commands run in the previous section, you should have obtained a ‘_resfinder’ output directory for each of the samples analysed. Out of the various output files, the following ones are key for interpretation:

- ```ResFinder_results_table.txt```: summary of detected acquired AMR genes by antibiotic class, including BLAST statistics such as percentage of nucleotide identity or percentage of gene length covered.
- ```PointFinder_table.txt```: summary of chromosomal genes scanned for AMR point mutations for the chosen bacterial species.
- ```PointFinder_results.txt```: detected AMR chromosomal point mutations in your sample, and associated phenotypic resistance.
- ```pheno_table_[species].txt```: ResFinder (including PointFinder) WGS-predicted phenotypes for the bacterial species chosen, including the detected genetic determinants supporting such prediction.
- ```pheno_table.txt```: ResFinder (including PointFinder) WGS-predicted phenotypes for all antibiotics included in this database. These results should be interpreted with caution, and the file pheno_table_[species].txt prioritised for reporting.

**IMPORTANT**: Out of all these output files, we should focus on ```pheno_table_[species].txt```, as it contains WGS-predicted phenotypes that are specific for the chosen bacterial species. Also, look at file PointFinder_results.txt to extract AMR point mutations.
