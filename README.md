# Scripts used for data analysis in the MARPLE _Pgt_ paper

## Variant Calling
For the analysis of genomic samples ([BWA](https://bio-bwa.sourceforge.net/) version 0.7.5; [Samtools](https://www.htslib.org/) version 0.1.19; [bamaddrg](https://github.com/ekg/bamaddrg) version 0.1.0):
```bash
bwa mem -M -t 8 [Reference.fna] [Sample_R1.fastq Sample_R2.fastq] | samtools view -S -b > [Sample.bam]
bamaddrg -b [Sample.bam] -s "sample" -r "sample" > [Sample_unsorted.bam]
samtools sort [Sample_unsorted.bam] [Sample_sorted.bam]
samtools index [Sample_sorted.bam]
samtools mpileup -f [Reference.fna] [Sample_sorted.bam] | gzip -9 -c > [Sample.pileup.gz]
```

For the analysis of transcriptomic samples ([STAR](https://github.com/alexdobin/STAR) version 2.5.a; [Samtools](https://www.htslib.org/) version 0.1.19; [GATK](https://gatk.broadinstitute.org/hc/en-us) version 4.0):
```bash
STAR --genomeDir [Reference/] --readFilesIn [Sample_R1.fastq Sample_R2.fastq] --outSAMtype BAM SortedByCoordinate
STAR --runMode genomeGenerate --genomeDir [Reference/] --genomeFastaFiles [Reference.fna] --sjdbGTFtagExonParentTranscript "Parent" --sjdbGTFfile [Reference.gff]
samtools index [Sample.bam]
gatk SplitNCigarReads -R [Reference.fna] -I [Sample.bam] -O [Sample_trimmed_aligned.bam]
samptools mpileup -f [Reference.fna] [Sample_trimmed_aligned.bam] | gzip -9 -c > [Sample.pileup.gz]
```

Obtain SNP information and Allele frequencies from the resulting pileup files:
```bash
gunzip -c [Sample.pileup.gz] | python compsnps_sampileup.py > [Sample_SNP_ratios.txt]
cat [Sample_SNP_ratios.txt] | perl compsnps_filter.pl 15 0.2 0.8 > [Sample_SNP_freq_15x.txt]
perl tab2ggplot2input_allele_freq.pl "sample" [Sample_SNP_freq_15x.txt] [OutDir/]
```

## _Pst_ Gene Amplification Checks
+ To get the percentage of genes that amplified in the _Pst_ samples using MARPLE, the aligned sequence file was analysed using a python script (`check-gene-coverage.py`).
+ After obtaining the percentage of gene amplification for the polymorphic genes, the data was visualised in a barchart using R (`GeneAmpl.R`).

## Select Polymorphic Genes
To perform the analysis on the SNP/b for each gene from the isolates included in the global population analysis, need to start from the consensus sequences:
```
python filter_and_sample_genes.py [Reference.gff] [alignments_OutDir/] [reports_OutDir/] [consensus/] [outgroup_consensus/]
```
Then, use [TreeCMP](https://github.com/TreeCmp/TreeCmp) to analyse the similarity between the main trees at each threshold and the 60% subsets (`treecmp_check_RF.sh`)

## Chromosome-level Gene Visualisation in Circos Plot
To map the polymorphic genes in _Pgt_ against the chromosome-level genome assembly available for the pathogen ([RKQQC](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002762355.2/)), we first run a BLAST search of the gene sequence to look for the corresponding genes in RKQQC, which were assigned to a chromosome.

The BLAST output is then analysed to look for genes with a sequence similarity >95% and bit-score >1,000, and for genes with >99% similarity that have another identical hit in a different chromosome.

The jupyter notebook to run the analysis: `ChromBLAST_CircosPlot.ipynb`

> [!NOTE]
> The script uses the `sampled_cds.csv` file, created during the gene selection step (_vide supra_), to plot the SNP/b percentage in a new panel in the Circos plot.


