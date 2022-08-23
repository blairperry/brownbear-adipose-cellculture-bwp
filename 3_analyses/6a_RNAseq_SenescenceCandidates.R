
library(tidyverse)
library(edgeR)
library(pheatmap)
library(colorspace)

gene_all <- read.csv("2_data/ell_culture_gene_count_matrix.csv",header=T, row.names=1)

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

normCounts <- cpm(y)

# Read in Swissprot Best Blast results
bblast <- read_csv('2_data/IDmap_2019_02_20.csv') %>%
  select(gene=1,human_id=3) %>%
  filter(!is.na(human_id)) %>%
  # mutate(human_id = str_split_fixed(human_id,'[_]',2)[,1]) %>%
  filter(str_detect(gene,'gene')) %>%
  left_join(id.convert,by=c('human_id'='id'))

normCounts.newIDs <- normCounts %>%
  as.data.frame() %>%
  rownames_to_column('old_id') %>%
  mutate(old_id = str_split_fixed(old_id,'[|]',2)[,1]) %>%
  left_join(bblast,by=c('old_id'='gene')) %>%
  filter(new_id != '<NA>') %>%
  select(-human_id) %>%
  select(-old_id) %>%
  select(new_id,everything())


# Read in senescence genes ------------------------------------------------
sen.candgenes <- read_tsv('2_data/SenescenceGeneList_06.15.22.txt')


# Get senescence gene expression table and format for heatmaps ------------
sen.normCounts <- normCounts.newIDs %>%
  filter(new_id %in% sen.candgenes$swissprot) %>%
  left_join(sen.candgenes,by=c('new_id'='swissprot')) %>%
  mutate(new_id = make.names(new_id, unique = TRUE),
         symbol = make.names(symbol,unique = TRUE))

sen.heatdata <- sen.normCounts %>%
  mutate(notes = factor(notes,levels=c('upregulated','sasp','downregulated'))) %>%
  arrange(notes,symbol) %>%
  select(-prot_id,-notes,-new_id) %>%
  select(symbol,everything()) %>%
  column_to_rownames('symbol')

sen.annot.row <- sen.normCounts %>%
  select(symbol,notes) %>%
  column_to_rownames('symbol')

sen.annot.col <- as.data.frame(names(sen.normCounts)) %>%
  select(sample=1) %>%
  filter(str_detect(sample,'symbol|prot|notes|new_id',negate = T)) %>%
  mutate(treatment = str_sub(sample,2)) %>%
  column_to_rownames('sample')


paletteLength <- 50
lim <-  3
myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(sen.heatdata,
         scale = 'row',
         cluster_rows = F,
         gaps_row = c(10,24),
         border_color = NA,
         color = diverging_hcl(palette = 'Blue Red 3',n = 50),
         breaks=myBreaks,
         cutree_cols = 3,
         annotation_row = sen.annot.row,
         annotation_col = sen.annot.col)


####


# Overlap DE results with senescence candidate genes ----------------------


de.HA_HH <- read_csv('2_data/de_gene_lists/tt_HA_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id) %>% filter(new_id %in% sen.candgenes$swissprot)
de.AA_HH <- read_csv('2_data/de_gene_lists/tt_AA_HH.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)
de.HPD_HH <- read_csv('2_data/de_gene_lists/tt_HPD_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)
de.AH_HH <- read_csv('2_data/de_gene_lists/tt_AH_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)
de.AA_AH <- read_csv('2_data/de_gene_lists/tt_AA_AH.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)
de.AH_APD <- read_csv('2_data/de_gene_lists/tt_AH_APD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot) %>% ## Flipping the sign of the FC here so that the Hib-serum treatment is second
  mutate(logFC = -1 * logFC) %>%
  mutate(temp = treat1,
         treat1 = treat2,
         treat2 = temp) %>%
  select(-temp)
de.HA_HPD <- read_csv('2_data/de_gene_lists/tt_HA_HPD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)
de.AA_APD <- read_csv('2_data/de_gene_lists/tt_AA_APD.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)%>% filter(new_id %in% sen.candgenes$swissprot)

all.deRes <- de.HA_HH %>%
  bind_rows(de.AA_HH ,
            de.HPD_HH,
            de.AH_HH ,
            de.AA_AH ,
            de.AH_APD,
            de.HA_HPD,
            de.AA_APD) %>%
  left_join(sen.candgenes,by=c('new_id'='swissprot')) %>%
  mutate(comparison = paste(treat1,' vs. ',treat2,sep = '')) %>%
  mutate(notes = str_replace_all(notes,'dowregulated','downregulated')) %>%
  mutate(notes = str_to_title(notes),
         notes = factor(notes,levels = c('Upregulated','Downregulated','Sasp')))

ggplot(all.deRes,aes(x=comparison,y=symbol,size=abs(logFC),fill=logFC)) +
  geom_point(pch=21) +
  scale_size_continuous(range = c(5,10)) +
  scale_fill_continuous_diverging(palette='Blue Red 3')+
  facet_grid(rows = vars(notes),scale='free_y',space = 'free') +
  theme_linedraw(base_size = 16) + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
