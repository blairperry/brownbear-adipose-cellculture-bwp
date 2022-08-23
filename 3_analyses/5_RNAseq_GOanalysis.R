
library(tidyverse)
library(WebGestaltR)
library(edgeR)


gene_all <- read_csv('2_data/brownbear-adipose-cellculture/Cell_culture_gene_count_matrix.csv')

# Filter BG set to be all genes used for pairwise comparisons

# Remove all annotated rRNA genes
gene_all = gene_all[-which(gene_all$gene_id %in% c("gene10974|LOC113259291","gene17084|LOC113266004","gene26915|LOC113246229","gene26916|LOC113246230","gene27720|LOC113247006","gene29258|LOC113248751","gene29259|LOC113248750","gene29260|LOC113248754","gene29262|LOC113248752","gene29263|LOC113248755","gene29856|LOC113249424","gene6633|LOC113254166","gene946|LOC113257075","gene947|LOC113257125","rna54509","rna54511")),]

# Trim out any transcripts with no counts
total <- gene_all[rowSums(gene_all[,-1]) > 0, ]

# Remove genes that are expressed in few individuals
keep <- rowSums(cpm(gene_all[,-1])>0.5) >= 3
length(which(keep == TRUE))

# Keep data that meet the filtering criteria
gene.counts <- gene_all[keep,]

# Read in uniprot-swissprot ID conversion table
id.convert <- read_tsv('2_data/brownbear-adipose-cellculture/uniprot_idConvertTable.txt') %>% 
  select(id=1,new_id=2)

# Read in Swissprot Best Blast results
bblast <- read_csv('2_data/brownbear-adipose-cellculture/IDmap_2019_02_20.csv') %>% 
  select(gene=1,human_id=3) %>%
  filter(!is.na(human_id)) %>% 
  # mutate(human_id = str_split_fixed(human_id,'[_]',2)[,1]) %>% 
  filter(str_detect(gene,'gene')) %>% 
  left_join(id.convert,by=c('human_id'='id'))


# Make background gene list for GO analyses

test.bg <- gene.counts %>% 
  select(gene_id) %>% 
  mutate(gene = str_split_fixed(gene_id,'[|]',2)[,1]) %>% 
  left_join(bblast,by=c('gene'='gene')) %>% 
  select(1,2,entry=3,human_id=4)

# write_tsv(test.bg[,'human_id'],'~/Downloads/testBG.txt')

# Write function to run GO analysis ---------

runGO <- function(infile,bgfile) {
  
    go.res.bp <- WebGestaltR(enrichMethod = 'ORA',
                               enrichDatabase = c('geneontology_Biological_Process_noRedundant'),
                               interestGene = infile$human_id,
                               interestGeneType = 'uniprotswissprot',
                               referenceGene = test.bg$human_id,
                               referenceGeneType = 'uniprotswissprot',
                               sigMethod = 'top',
                               minNum=5,
                               maxNum = 2000,
                               topThr = 500,
                               isOutput = F) %>% 
      mutate(database = 'biological_process')
    go.res.cc <- WebGestaltR(enrichMethod = 'ORA',
                                  enrichDatabase = c('geneontology_Cellular_Component_noRedundant'),
                                  interestGene = infile$human_id,
                                  interestGeneType = 'uniprotswissprot',
                                  referenceGene = test.bg$human_id,
                                  referenceGeneType = 'uniprotswissprot',
                                  sigMethod = 'top',
                                  minNum=5,
                                  maxNum = 2000,
                                  topThr = 500,
                                  isOutput = F)%>% 
      mutate(database = 'cellular_component')
    go.res.mf <- WebGestaltR(enrichMethod = 'ORA',
                                  enrichDatabase = c('geneontology_Molecular_Function_noRedundant'),
                                  interestGene = infile$human_id,
                                  interestGeneType = 'uniprotswissprot',
                                  referenceGene = test.bg$human_id,
                                  referenceGeneType = 'uniprotswissprot',
                                  sigMethod = 'top',
                                  minNum=5,
                                  maxNum = 2000,
                                  topThr = 500,
                                  fdrThr = 0.05,
                                  isOutput = F)%>% 
      mutate(database = 'molecular_function')
  
    go.res <- go.res.bp %>% 
      bind_rows(go.res.cc,go.res.mf) %>% 
      filter(FDR < 0.05)
  
    return(go.res)
}

# Read in all pairwise comparison analyses --------------------------------

de.HA_HH <- read_csv('2_data/de_gene_lists/tt_HA_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.AA_HH <- read_csv('2_data/de_gene_lists/tt_AA_HH.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.HPD_HH <- read_csv('2_data/de_gene_lists/tt_HPD_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.AH_HH <- read_csv('2_data/de_gene_lists/tt_AH_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.AA_AH <- read_csv('2_data/de_gene_lists/tt_AA_AH.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.AH_APD <- read_csv('2_data/de_gene_lists/tt_AH_APD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.HA_HPD <- read_csv('2_data/de_gene_lists/tt_HA_HPD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
de.AA_APD <- read_csv('2_data/de_gene_lists/tt_AA_APD.csv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)

# de.AA_HPD <- read_csv('2_data/de_gene_lists/tt_AA_HPD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
# de.AH_HA <- read_csv('2_data/de_gene_lists/tt_AH_HA.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
# de.AH_HPD <- read_csv('2_data/de_gene_lists/tt_AH_HPD.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
# de.APD_HA <- read_csv('2_data/de_gene_lists/tt_APD_HA.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)
# de.APD_HH <- read_csv('2_data/de_gene_lists/tt_APD_HH.csv')  %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id)

write_tsv(de.AH_APD,'3_analyses/forShinyGOTest_AAHsHG.tsv')


# Run all GO analyses -----------------------------------------------------

goRes.de.HA_HH  <- runGO(de.HA_HH  ,test.bg) %>% mutate(comparison = "HA_HH")
goRes.de.AA_HH  <- runGO(de.AA_HH  ,test.bg) %>% mutate(comparison = "AA_HH")
goRes.de.HPD_HH <- runGO(de.HPD_HH ,test.bg) %>% mutate(comparison = "HPD_HH")
goRes.de.AH_HH  <- runGO(de.AH_HH  ,test.bg) %>% mutate(comparison = "AH_HH")
goRes.de.AH_APD <- runGO(de.AH_APD ,test.bg) %>% mutate(comparison = "AH_APD")
goRes.de.HA_HPD <- runGO(de.HA_HPD ,test.bg) %>% mutate(comparison = "HA_HPD")
goRes.de.AA_APD <- runGO(de.AA_APD ,test.bg) %>% mutate(comparison = "AA_APD") 
goRes.de.AA_AH  <- runGO(de.AA_AH  ,test.bg) %>% mutate(comparison = "AA_AH")


goRes.all <- goRes.de.HA_HH %>% 
  bind_rows(goRes.de.AA_HH,goRes.de.HPD_HH,goRes.de.AH_HH,goRes.de.AH_APD,goRes.de.HA_HPD,goRes.de.AA_APD,goRes.de.AA_AH) %>% 
  mutate(comparison = str_replace_all(comparison,'_',' vs. '))

goRes.all %>% 
  group_by(comparison) %>% 
  tally()

# write_csv(goRes.all,'3_analyses/__NewGOResultSuppTable_06.13.22.csv')

goRes.focal <- goRes.de.HA_HH %>% 
  bind_rows(goRes.de.AA_HH,goRes.de.HPD_HH)

# Plot overlap in focal comparison GO results
ggplot(goRes.focal,aes(x=comparison,y=description,size=enrichmentRatio)) +
  geom_point()



# Summarize and compare with Rrvigo ----------------------------------------

library(rrvgo)

simMatrix.bp <- calculateSimMatrix(goRes.focal$geneSet,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms.bp <- reduceSimMatrix(simMatrix.bp,
                                threshold=0.6,
                                orgdb="org.Hs.eg.db")
reducedTerms.bp.simple <- reducedTerms.bp %>% select(go,parentTerm)

simMatrix.cc <- calculateSimMatrix(goRes.focal$geneSet,
                                   orgdb="org.Hs.eg.db",
                                   ont="CC",
                                   method="Rel")

reducedTerms.cc <- reduceSimMatrix(simMatrix.cc,
                                   threshold=0.6,
                                   orgdb="org.Hs.eg.db")
reducedTerms.cc.simple <- reducedTerms.cc %>% select(go,parentTerm)

simMatrix.mf <- calculateSimMatrix(goRes.focal$geneSet,
                                   orgdb="org.Hs.eg.db",
                                   ont="MF",
                                   method="Rel")

reducedTerms.mf <- reduceSimMatrix(simMatrix.mf,
                                   threshold=0.6,
                                   orgdb="org.Hs.eg.db")
reducedTerms.mf.simple <- reducedTerms.mf%>% select(go,parentTerm)

reducedTerms.all.simple <- reducedTerms.bp.simple %>% bind_rows(reducedTerms.cc.simple) %>% unique()

goRes.all.parent <- goRes.focal %>%
  left_join(reducedTerms.all.simple,by = c('geneSet'='go'))



##

# Plot comparison of GO terms ---------------------------------------------
library(pheatmap)

ggplot(goRes.all.parent %>% filter(!is.na(parentTerm)),aes(x=comparison,y=parentTerm,fill=database)) +
  geom_point(color='black',pch=21,size=6) +
  facet_grid(rows=vars(database),space = 'free_y',scale='free_y') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

# goRes.bp.heatdata <- goRes.all.parent %>%
#   filter(!is.na(parentTerm)) %>%
#   filter(str_detect(database,'rocess')) %>%
#   select(comparison,parentTerm) %>%
#   unique() %>%
#   mutate(val = 1) %>%
#   pivot_wider(names_from = comparison,values_from = val,values_fill = 0) %>%
#   as.data.frame() %>%
#   column_to_rownames('parentTerm')
# 
# pheatmap(goRes.bp.heatdata,
#          border_color = 'black',
#          treeheight_row = 10,treeheight_col = 10,
#          cellheight = 10,cellwidth = 10,
#          color = c('aliceblue','dodgerblue4')
# )


###


# WGCNA Module GO Analysis ------------------------------------------------

wgcna.mods <- read_tsv('3_analyses/WGCNA_GeneModuleMembership_05.26.22.tsv') %>% mutate(gene = str_split_fixed(gene,'[|]',2)[,1]) %>% left_join(bblast,by=c('gene'='gene')) %>% mutate(human_id = new_id) %>% select(-new_id)

# Running GO analyses for all modules that have at least one module-trait correlation > 0.5 with p-val < 0.05
purple.go.res <- runGO(wgcna.mods %>% filter(module_color=='purple')  ,wgcna.mods)  %>% mutate(module="purple")
green.go.res <- runGO(wgcna.mods %>% filter(module_color=='green')  ,wgcna.mods) %>% mutate(module="green")
red.go.res <- runGO(wgcna.mods %>% filter(module_color=='red')  ,wgcna.mods) %>% mutate(module="red")
black.go.res <- runGO(wgcna.mods %>% filter(module_color=='black')  ,wgcna.mods) %>% mutate(module="black")
pink.go.res <- runGO(wgcna.mods %>% filter(module_color=='pink')  ,wgcna.mods) %>% mutate(module="pink")
midnightblue.go.res <- runGO(wgcna.mods %>% filter(module_color=='midnightblue')  ,wgcna.mods) %>% mutate(module="midnightblue")
darkgrey.go.res <- runGO(wgcna.mods %>% filter(module_color=='darkgrey')  ,wgcna.mods) %>% mutate(module="darkgrey")
darkturquoise.go.res <- runGO(wgcna.mods %>% filter(module_color=='darkturquoise')  ,wgcna.mods) %>% mutate(module="darkturquoise")
brown.go.res <- runGO(wgcna.mods %>% filter(module_color=='brown')  ,wgcna.mods) %>% mutate(module="brown")
darkred.go.res <- runGO(wgcna.mods %>% filter(module_color=='darkred')  ,wgcna.mods) %>% mutate(module="darkred")
cyan.go.res <- runGO(wgcna.mods %>% filter(module_color=='cyan')  ,wgcna.mods) %>% mutate(module="cyan")
lightcyan.go.res <- runGO(wgcna.mods %>% filter(module_color=='lightcyan')  ,wgcna.mods) %>% mutate(module="lightcyan")
navajowhite.go.res  <- runGO(wgcna.mods %>% filter(module_color=='navajowhite2')  ,wgcna.mods) %>% mutate(module="navajowhite")
blue.go.res  <- runGO(wgcna.mods %>% filter(module_color=='blue')  ,wgcna.mods) %>% mutate(module="blue")

all.wgcna.go.res <- purple.go.res %>% 
  bind_rows(green.go.res, red.go.res, black.go.res, pink.go.res, midnightblue.go.res, darkgrey.go.res, darkturquoise.go.res, brown.go.res, darkred.go.res, cyan.go.res, lightcyan.go.res, navajowhite.go.res, blue.go.res)



wgcna.simMatrix.bp <- calculateSimMatrix(all.wgcna.go.res$geneSet,
                                   orgdb="org.Hs.eg.db",
                                   ont="BP",
                                   method="Rel")

wgcna.reducedTerms.bp <- reduceSimMatrix(wgcna.simMatrix.bp,
                                   threshold=0.6,
                                   orgdb="org.Hs.eg.db")
wgcna.reducedTerms.bp.simple <- wgcna.reducedTerms.bp %>% select(go,parentTerm)

wgcna.simMatrix.cc <- calculateSimMatrix(all.wgcna.go.res$geneSet,
                                   orgdb="org.Hs.eg.db",
                                   ont="CC",
                                   method="Rel")

wgcna.reducedTerms.cc <- reduceSimMatrix(wgcna.simMatrix.cc,
                                   threshold=0.6,
                                   orgdb="org.Hs.eg.db")
wgcna.reducedTerms.cc.simple <- wgcna.reducedTerms.cc %>% select(go,parentTerm)

wgcna.simMatrix.mf <- calculateSimMatrix(all.wgcna.go.res$geneSet,
                                   orgdb="org.Hs.eg.db",
                                   ont="MF",
                                   method="Rel")

wgcna.reducedTerms.mf <- reduceSimMatrix(wgcna.simMatrix.mf,
                                   threshold=0.6,
                                   orgdb="org.Hs.eg.db")

wgcna.reducedTerms.mf.simple <- wgcna.reducedTerms.mf%>% select(go,parentTerm)

wgcna.reducedTerms.all.simple <- wgcna.reducedTerms.mf %>% bind_rows(wgcna.reducedTerms.cc.simple,wgcna.reducedTerms.bp.simple) %>% unique()

wgcna.goRes.all.parent <- all.wgcna.go.res %>%
  left_join(wgcna.reducedTerms.all.simple,by = c('geneSet'='go'))

ggplot(wgcna.goRes.all.parent %>% filter(parentTerm != 'NA'),aes(x=module,y=parentTerm)) +
  geom_point() +
  facet_grid(rows=vars(database),scales = 'free_y',space = 'free')

wgcna.goRes.bp.parent.heatdata <- wgcna.goRes.all.parent %>% 
    filter(!is.na(parentTerm)) %>%
    filter(str_detect(database,'process')) %>% 
    select(module,parentTerm) %>%
    unique() %>% 
    mutate(val = 1) %>%
    pivot_wider(names_from = module,values_from = val,values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames('parentTerm')

wgcna.goRes.cc.parent.heatdata <- wgcna.goRes.all.parent %>% 
  filter(!is.na(parentTerm)) %>%
  filter(str_detect(database,'cellular')) %>% 
  select(module,parentTerm) %>%
  unique() %>% 
  mutate(val = 1) %>%
  pivot_wider(names_from = module,values_from = val,values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames('parentTerm')

wgcna.goRes.mf.parent.heatdata <- wgcna.goRes.all.parent %>% 
  filter(!is.na(parentTerm)) %>%
  filter(str_detect(database,'molec')) %>% 
  select(module,parentTerm) %>%
  unique() %>% 
  mutate(val = 1) %>%
  pivot_wider(names_from = module,values_from = val,values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames('parentTerm')

pheatmap::pheatmap(wgcna.goRes.bp.parent.heatdata,
         border_color = 'black',
         treeheight_row = 10,treeheight_col = 10,
         cellheight = 10,cellwidth = 10,
         main='Biological Process\n',
         color = c('white','dodgerblue4')
)

pheatmap::pheatmap(wgcna.goRes.cc.parent.heatdata,
                   border_color = 'black',
                   treeheight_row = 10,treeheight_col = 10,
                   cellheight = 10,cellwidth = 10,
                   main='Cellular Component\n',
                   color = c('white','dodgerblue4')
)

pheatmap::pheatmap(wgcna.goRes.mf.parent.heatdata,
                   border_color = 'black',
                   treeheight_row = 10,treeheight_col = 10,
                   cellheight = 10,cellwidth = 10,
                   main='Molecular Function\n',
                   color = c('white','dodgerblue4')
)







