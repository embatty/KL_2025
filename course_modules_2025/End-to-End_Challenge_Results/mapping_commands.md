
```
cd /home/data
mkdir mapping
```

```
unzip Fastq\ files\ for\ End-to-End\ challenge\ samples-20250121T100355Z-001.zip
cp Fastq\ files\ for\ End-to-End\ challenge\ samples/cpe107_* ./mapping/
cp Fastq\ files\ for\ End-to-End\ challenge\ samples/cpe108_* ./mapping/
```

```
cd mapping
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR538/002/ERR5386322/ERR5386322_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR538/002/ERR5386322/ERR5386322_2.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR538/008/ERR5386328/ERR5386328_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR538/008/ERR5386328/ERR5386328_2.fastq.gz
```

```
fastp -i ERR5386322_1.fastq.gz -I ERR5386322_2.fastq.gz -o out.ERR5386322_1.fastq.gz -O out.ERR5386322_2.fastq.gz
fastp -i ERR5386328_1.fastq.gz -I ERR5386328_2.fastq.gz -o out.ERR5386328_1.fastq.gz -O out.ERR5386328_2.fastq.gz
fastp -i cpe107_R1.fastq.gz -I cpe107_R2.fastq.gz -o out.cpe107_R1.fastq.gz -O out.cpe107_R2.fastq.gz
fastp -i cpe108_R1.fastq.gz -I cpe108_R2.fastq.gz -o out.cpe108_R1.fastq.gz -O out.cpe108_R2.fastq.gz
```

```
conda deactivate
conda activate snippy
snippy --cpus 4 --outdir ERR5386322_snippy --reference cpe070_Kpn-ST11-NDM1.chr.fasta --R1 out.ERR5386322_1.fastq.gz --R2 out.ERR5386322_2.fastq.gz
snippy --cpus 4 --outdir ERR5386328_snippy --reference cpe070_Kpn-ST11-NDM1.chr.fasta --R1 out.ERR5386328_1.fastq.gz --R2 out.ERR5386328_2.fastq.gz
snippy --cpus 4 --outdir cpe107_snippy --reference cpe070_Kpn-ST11-NDM1.chr.fasta --R1 out.cpe107_R1.fastq.gz --R2 out.cpe107_R2.fastq.gz
snippy --cpus 4 --outdir cpe108_snippy --reference cpe070_Kpn-ST11-NDM1.chr.fasta --R1 out.cpe108_R1.fastq.gz --R2 out.cpe108_R2.fastq.gz
conda deactivate
conda activate amr
```



