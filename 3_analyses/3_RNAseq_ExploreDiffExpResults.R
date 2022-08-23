
library(tidyverse)
library(eulerr)
library(viridis)
library(edgeR)
library(pheatmap)
library(colorspace)



gene_all <- read.csv("./2_data/Cell_culture_gene_count_matrix.csv",header=T, row.names=1)

# Remove all annotated rRNA genes
gene_all = gene_all[-which(rownames(gene_all) %in% c("gene10974|LOC113259291","gene17084|LOC113266004","gene26915|LOC113246229","gene26916|LOC113246230","gene27720|LOC113247006","gene29258|LOC113248751","gene29259|LOC113248750","gene29260|LOC113248754","gene29262|LOC113248752","gene29263|LOC113248755","gene29856|LOC113249424","gene6633|LOC113254166","gene946|LOC113257075","gene947|LOC113257125","rna54509","rna54511")),]


# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all) > 0, ]

# Remove genes that are expressed in few individuals
keep <- rowSums(cpm(gene_all)>0.5) >= 3
length(which(keep == TRUE))

# Keep data that meet the filtering criteria
gene.counts <- gene_all[keep,]

#Design Matrices with bear as a blocking factor
group = factor(c("AA", "AH", "AG", "HA", "HH", "HG","AA", "AH", "AG", "HA", "HH", "HG","AA", "AH", "AG", "HA", "HH", "HG","AA", "AH", "AG", "HA", "HH", "HG","AA", "AH", "AG", "HA", "HH", "HG","AA", "AG", "HA", "HH", "HG"))
y <- DGEList(counts=gene.counts,group=group)
y <- calcNormFactors(y)

Bears.dex <- factor(c("Cooke","Cooke","Cooke","Cooke","Cooke","Cooke","Frank","Frank","Frank","Frank","Frank","Frank","John","John","John","John","John","John","Oakley","Oakley","Oakley","Oakley","Oakley","Oakley","Pacino","Pacino","Pacino","Pacino","Pacino","Pacino","Roan","Roan","Roan","Roan","Roan"))
Cell.Serum <- factor(c("A.A","A.H","A.G","H.A","H.H","H.G","A.A","A.H","A.G","H.A","H.H","H.G","A.A","A.H","A.G","H.A","H.H","H.G","A.A","A.H","A.G","H.A","H.H","H.G","A.A","A.H","A.G","H.A","H.H","H.G","A.A","A.G","H.A","H.H","H.G"))
fat.design <- model.matrix(~ 0 + Cell.Serum + Bears.dex, data = y$samples)

#Estimate Dispersions
y <- estimateDisp(y, fat.design, robust=TRUE)
y$common.dispersion


# Read in result files  -----------------

tt_AA.AH <- read_csv('2_data/de_gene_lists/tt_AA_AH.csv') 
tt_AA.APD <- read_csv('2_data/de_gene_lists/tt_AA_APD.csv') 
tt_AA.HH <- read_csv('2_data/de_gene_lists/tt_AA_HH.csv') 
tt_AA.HPD <- read_csv('2_data/de_gene_lists/tt_AA_HPD.csv')  
tt_AH.APD <- read_csv('2_data/de_gene_lists/tt_AH_APD.csv')  
tt_AH.HA <- read_csv('2_data/de_gene_lists/tt_AH_HA.csv')  
tt_AH.HH <- read_csv('2_data/de_gene_lists/tt_AH_HH.csv')  
tt_AH.HPD <- read_csv('2_data/de_gene_lists/tt_AH_HPD.csv')  
tt_APD.HA <- read_csv('2_data/de_gene_lists/tt_APD_HA.csv')  
tt_APD.HH <- read_csv('2_data/de_gene_lists/tt_APD_HH.csv')  
tt_HA.HH <- read_csv('2_data/de_gene_lists/tt_HA_HH.csv')  
tt_HA.HPD <- read_csv('2_data/de_gene_lists/tt_HA_HPD.csv')  
tt_HPD.HH <- read_csv('2_data/de_gene_lists/tt_HPD_HH.csv')  


# Combine and summarize all DE --------------------------------------------

tt_all <- tt_HPD.HH %>% 
  bind_rows(tt_HA.HPD,tt_HA.HH,tt_APD.HH,tt_APD.HA,tt_AH.HPD,tt_AH.HH,tt_AH.HA,tt_AH.APD,tt_AA.HPD,tt_AA.HH,tt_AA.APD,tt_AA.AH)

tt_all.summary <- tt_all %>% 
  group_by(treat1,treat2) %>% 
  tally() %>% 
  ungroup() %>% 
  complete(treat1,treat2) 

tt_all.summary.upDown <- tt_all %>% 
  group_by(direction,treat1,treat2) %>% 
  tally() %>% 
  ungroup(treat1,treat2) %>% 
  complete(treat1,treat2) %>% 
  mutate(direction = factor(direction,levels=c('Upregulated','Downregulated')),
         n = ifelse(is.na(n),0,n))

tt_all.summary.upDown %>% 
  filter(n>0) %>% 
  ggplot(aes(x=treat2,y=treat1,size=n,fill=direction)) +
  geom_point(pch=21,color='black',alpha=0.6) +
  geom_text(data= . %>% filter(n>100),aes(label=n),color='black',size=2.5) +
  geom_text(data= . %>% filter(n<100),aes(label=n),color='black',size=2.5,nudge_x = 0.15,nudge_y=0.2) +
  scale_size_continuous(range = c(1,11)) +
  scale_fill_manual(values = c('Upregulated'='firebrick','Downregulated'='Dodgerblue')) +
  labs(x='Treatment 2',y='Treatment 1') +
  facet_wrap(~direction,ncol=1) +
  theme_linedraw()

tt_all.summary.upDown %>% 
  filter(n>0) %>%
  group_by(treat1,treat2) %>% 
  mutate(total = sum(n)) %>% 
  ggplot(aes(x=treat2,y=treat1,size=total)) +
  geom_point(pch=21,color='black') +
  geom_text(data= . %>% filter(n>100),aes(label=total),color='black',size=2.5) +
  geom_text(data= . %>% filter(n<100),aes(label=total),color='black',size=2.5,nudge_x = 0.15,nudge_y=0.2) +
  scale_size_continuous(range = c(1,11)) +
  labs(x='Treatment 2',y='Treatment 1') +
  # facet_wrap(~direction,ncol=1) +
  theme_linedraw()

#
    

# UPDATE 06.25.22 - Combine and summarize all FOCAL comparisons --------------------------------------------
library(eulerr)
set.seed(1)

comp.levels <- c('AA vs. HH','AH vs. HH','AA vs. HA','AG vs. HG','HA vs. HH','HG vs. HH','AA vs. AG','HA vs. HG','AA vs. AH','AH vs. AG')

tt_focal <- tt_AH.HH %>% 
  bind_rows(tt_AA.HH,tt_HA.HH,tt_HPD.HH,tt_AA.APD,tt_AA.AH,tt_AH.APD,tt_HA.HPD) %>% 
  mutate(comparison = paste(treat1, ' vs. ',treat2,sep = '')) %>% 
  mutate(comparison = factor(comparison, levels=comp.levels))


tt_focal.summary <- tt_focal %>% 
  group_by(comparison,direction) %>% 
  tally() %>% 
  ungroup() %>% 
  add_row(comparison = 'AA vs. HA',direction = 'Upregulated',n=0) %>%  # Manually add entries for comparisons with zero DE genes
  add_row(comparison = 'AA vs. HA',direction = 'Downregulated',n=0) %>% 
  add_row(comparison = 'AG vs. HG',direction = 'Upregulated',n=0) %>% 
  add_row(comparison = 'AG vs. HG',direction = 'Downregulated',n=0) %>% 
  mutate(direction = factor(direction,levels=c('Upregulated','Downregulated'))) %>% 
  group_by(comparison) %>% 
  mutate(total = sum(n))

ggplot(tt_focal.summary,aes(x=n,y=reorder(comparison,total))) +
  geom_bar(stat='identity',orientation = 'y',fill='grey20') +
  geom_text(data=tt_focal.summary %>% filter(direction=='Upregulated'),aes(label=total,x=total+500),size=4) +
  labs(x='# of DE Genes',y='Comparison',fill='Direction:') + 
  scale_x_continuous(expand=c(0,0),limits=c(0,7500)) +
  theme_classic(base_size = 16) + theme(legend.position = 'bottom')

# Euler venns of key overlaps

vennData.all <- list(`HA vs. HH` = tt_HA.HH$gene, 
                     `HG vs. HH` = tt_HPD.HH$gene, 
                     `AA vs. HH` = tt_AA.HH$gene, 
                     `AH vs. HH` = tt_AH.HH$gene,
                     `AA vs. AH` = tt_AA.AH$gene,
                     `AH vs. AG` = tt_AH.APD$gene,
                     `HA vs. HG` = tt_HA.HPD$gene,
                     `AA vs. AG` = tt_AA.APD$gene,
                     `AG vs. HG` = c(),
                     `AA vs. HA` = c())

plot(euler(vennData.all[c('HA vs. HH','HG vs. HH')]), quantities = TRUE,alpha=0.7,edges='black',labels=T) # A and G serum elicit similar response in hib cells
plot(euler(vennData.all[c('HA vs. HH','HG vs. HH','AA vs. HH')]), quantities = TRUE,alpha=0.5,edges='black',labels=T,fill=magma(3)) # Large overlap in all "active-ish" treatments compared to HH
plot(euler(vennData.all[c('AA vs. HH','HA vs. HH')]), quantities = TRUE,alpha=0.7,edges='black',labels=T) # AA and HA have large overlap in DE genes compared to HH
plot(euler(vennData.all[c('AH vs. HH','AA vs. HH')]), quantities = TRUE,alpha=0.7,edges='black',labels=T) # AH and AA have some overlap compared to HH

three.way.intersect <- as.data.frame(intersect(intersect(tt_HA.HH$gene,tt_AA.HH$gene),tt_HPD.HH$gene)) %>% 
  select(gene_id=1)


# Heatmap of genes from above venn diagram --------------------------------

normCounts <- cpm(y) %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  filter(gene %in% three.way.intersect$gene_id)
  

normCounts.heatdata <- normCounts %>%
  select(gene,contains('AA'),contains('HA'),contains('HG'),contains('HH')) %>%
  column_to_rownames('gene')

colAnnot <- as.data.frame(names(normCounts.heatdata)) %>% 
  select(sample=1) %>% 
  mutate(treatment = str_sub(sample,2)) %>% 
  column_to_rownames('sample')

pheatmap(log10(normCounts.heatdata+1),
         show_rownames = F,
         # scale = 'row',
         # cluster_cols = F,
         # cluster_rows = F,
         annotation_col = colAnnot,
         color = viridis(30),
         treeheight_row = 0,
         # gaps_col = c(6,12,18),
         cutree_cols = 2,
         border_color = NA
         )


# Assess overlap with DE genes identified in Jansen et al.  2021 ----------

cellCulture_HHvsAA <- read_csv('2_data/de_gene_lists/tt_AA_HH.csv') %>% 
  mutate(gene_id = str_split_fixed(gene,'[|]',2)[,1]) %>% 
  mutate(cc_direction = direction)

# Supplementary table 4 from Jansen et al. 2021
jansen_deGenes <- readxl::read_xlsx('../42003_2019_574_MOESM4_ESM.xlsx',skip = 1,na = 'NA') %>% 
  janitor::clean_names() %>% 
  filter(adipose_fdr < 0.05) %>% 
  mutate(t_direction = ifelse(adipose_log_fc > 0, 'Up','Down'))

overlap <- cellCulture_HHvsAA %>% filter(gene_id %in% jansen_deGenes$gene_id) %>% 
  left_join(jansen_deGenes,by = 'gene_id') %>% 
  select(gene, gene_id, logFC, adipose_log_fc, cc_direction, t_direction) %>% 
  mutate(symbol = str_split_fixed(gene, '[|]',2)[,2])

overlap %>% 
  group_by(cc_direction,t_direction) %>% 
  tally()

nrow(overlap) / nrow(cellCulture_HHvsAA) 

# 2225 (36.4%) DE genes were identified previously in bulk tissue RNAseq analyses 
# 659 are down in both, 736 are up in both
# 208 and 422 are down in cell and up in tissue, or the reverse, resepctively 

overlap.consistent <- overlap %>% 
  filter(cc_direction == t_direction)

# write_tsv(overlap.consistent,'../../3_analyses/OverlapDE_ConsistDir_05.09.22.tsv')


