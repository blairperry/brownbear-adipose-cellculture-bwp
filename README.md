# brownbear-adipose-cellculture-bwp
This repo contains code for data processing and analysis used in our iScience publication of experiments involving brown bear adipose cell culture stimulated with blood serum from active and hibernating bears. 

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
### Newly generated data:
- RNA-seq from from bear adipocyte cell culture
	- Raw data are available at NCBI Bioproject [PRJNA578991](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA578991)
- Proteomic data from bear blood serum
	- Raw and processed proteomic data are available at Panorama Public
		- ProteomeXchange ID: PXD023555
		- Access URL: [https://panoramaweb.org/grizzlybear.url](https://panoramaweb.org/grizzlybear.url)

## Contents
1. [Quality trimming, mapping, and processing of RNA-seq data](#1-brownbear-adipose-cellculture-bwp#1-quality-trimming-mapping-and-processing-of-rna-seq-data)
2. [Gene Expression Analyses - pairwise comparisons and related plotting](#2-gene-expression-analyses---pairwise-comparisons-and-related-plotting)
3. [Weighted Gene Correlation Network Analysis (WGCNA)](#3-weighted-gene-correlation-network-analysis-wgcna)
4. [Gene Ontology (GO) term overrepresentation analysis](#4-gene-ontology-go-term-overrepresentation-analysis)
5. [Candidate gene investigation](#5-candidate-gene-investigation)

---
### 1. Quality trimming, mapping, and processing of RNA-seq data
Note: These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak. Full SLURM scripts and additional detail can be found in [1_RNAseq_TrimAndMapping.rtf](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/1_RNAseq_TrimAndMapping.rtf).
- Assemblies used:
	- Brown bear: [NCBI Assembly ASM358476v1](https://www.ncbi.nlm.nih.gov/assembly/GCF_003584765.1/)

#### Quality assessment with FastQC and trimming with TrimGalore
```bash
# Run Fastqc
fastqc -o 1_FastQC/Results --noextract --nogroup bear_cell_culture/data/*fastq.gz

# Run TrimGalore
trim_galore --paired -q 20 --fastqc --fastqc_args "--noextract --nogroup --outdir 2_TrimGalore/fastqc/" --stringency 5 --illumina --length 50 -o trimmed_reads/ --clip_R1 12 --clip_R2 12 [path/to/read1] [path/to/read2]
```

#### Mapping using HiSat2
```bash
# Run HiSat2
hisat2 --phred33 -k 10 --met-file stats/${file_name}.stats --rg-id ${file_name} --rg SM:${file_name} --rg PL:illumina -p 1 --rna-strandness RF --fr --dta --un-conc-gz unmapped/${file_name}.unmapped -x Uarctoshorribilis/Uarctoshorribilis-genomic -1 [path/to/trimmed/read1] -2 [path/to/trimmed/read2] -S ${file_name}.sam

```

#### Fix mate information with Picard and sort SAM files with samtools
```bash
picard FixMateInformation INPUT=[path/to/samfile] OUTPUT=${file_name}.fix.sam VALIDATION_STRINGENCY=SILENT

samtools view -@ 10 -b -h ${file_name}.fix.sam > ${file_name}.bam

samtools sort -m 4G -o results/${file_name}_sorted_fixed.bam -O bam -T ${file_name} -@ 10 ${file_name}.bam
```

#### Collect mapping stats with Picard
```bash
picard CollectAlignmentSummaryMetrics R=Uarctoshorribilis/GCF_003584765.1_ASM358476v1_genomic.fna I=[path/to/sorted/bam] O=4_Fix_sort_bam/picard_mapping_stats/picard_output/${file_name}_alignmentsummarymetrics
```

#### Merge multiple bam files per sample with Picard
Each sample was sequenced with three runs, so need to merge the resulting bam files.
```bash
picard MergeSamFiles
I=[path/to/run1/sorted/bam] I=[path/to/run2/sorted/bam] I=[path/to/run3/sorted/bam] O=5_Merge_bams/results/${file_name}_temp.bam
  

samtools sort -m 4G -o 5_Merge_bams/results/${file_name}_merged.bam -O bam -T ${file_name} -@ 10 5_Merge_bams/results/${file_name}_temp.bam

rm 5_Merge_bams/results/*_temp.bam
```

#### Assemble transcripts and generate gene count matrix with Stringtie
```bash
# Assembly transcripts
stringtie [path/to/merged/bam] -o results/${file_name}.gtf -m 50 --rf -e -B -c 10 -p 2 Uarctoshorribilis/GCF_003584765.1_ASM358476v1_modified_mito.gff

# Generate gene count matrix with prepDE.py from Stringtie
python prepDE.py -i datatable -g 7_PrepDE/matrix_with_fixed_brown_bear/gene_count_matrix.csv -t 7_PrepDE/matrix_with_fixed_brown_bear/transcript_count_matrix.csv -l 88
# NOTE: `datatable` in above command has sample id in first column and path to the corresponding Stringtie result GTF file in the second column 

```

### 2. Gene Expression Analyses - Pairwise comparisons and related plotting
The following R script contains code used to:
- Normalize gene expression counts
- Perform differential expression analyses between treatment pairs
- Generate multi-dimensional scaling (MDS) plots to assess large-scale patterns in expression across samples
- Generate supplementary upset plot figures to compare differentially expressed genes across analyses

Link to Rscript: [2_RNAseq_DiffExpression.Rmd](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/2_RNAseq_DiffExpression.Rmd)

The following R script contains code used for subsequent post-hoc analysis and plotting of pairwise comparison results, inlcuding:
- Plotting total number of differentially expressed (DE) genes for focal analyses
- Plotting euler venn diagrams and heatmaps of overlapping DE genes between focal analyses
- Comparing DE genes with previous tissue level DE genes from Jansen et al. 2021

Link to Rscript: [3_RNAseq_ExploreDiffExpResults.R](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/3_RNAseq_ExploreDiffExpResults.R)

### 3. Weighted Gene Correlation Network Analysis (WGCNA)
The following R script contains code used to run WGCNA analyses and plot results

Link to Rscript: [4_RNAseq_WGCNA.R](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/tree/main/3_analyses#:~:text=5%20minutes%20ago-,4_RNAseq_WGCNA.R,-initial%20commit)

### 4. Gene Ontology (GO) term overrepresentation analysis
The following R script contains code used to run GO analyses on sets of DE genes and WGCNA modules

Link to Rscript: [5_RNAseq_GOanalysis.R](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/5_RNAseq_GOanalysis.R)

### 5. Candidate gene investigation
We curated lists of candidate genes associated with a) senescence phenotypes in adipocytes (based on Figure 5b in [Li et al. 2021](https://www.nature.com/articles/s41591-021-01501-8)) and b) the Complement and [Coagulation Cascades KEGG Pathway (hsa04610)](https://www.genome.jp/pathway/hsa04610) and investigate the expression of these candidate genes in our dataset. The following two R scripts contain code used to run these analyses.

Link to senescence Rscript: [6a_RNAseq_SenescenceCandidates.R](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/6a_RNAseq_SenescenceCandidates.R)

Link to complement Rscript: [6b_RNAseq_ComplementCandidates.R](https://github.com/blairperry/brownbear-adipose-cellculture-bwp/blob/main/3_analyses/6b_RNAseq_ComplementCandidates.R)
