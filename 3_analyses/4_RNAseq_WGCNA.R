
library(tidyverse)
library(WGCNA)
library(DESeq2)

# Import raw counts and normalize with VST  -------------------------------

rawCounts <- read_csv('2_data/Cell_culture_gene_count_matrix.csv') %>% 
  as.data.frame() %>% 
  column_to_rownames('gene_id')

sampleInfo <- as.data.frame(names(rawCounts)) %>% select(sample = 1) %>% 
  mutate(treatment = str_sub(sample,2))

# Set up DEseq data from count matrix, use drainage as design 
colData <- DataFrame(treatment= factor(sampleInfo$treatment))

dds <- DESeqDataSetFromMatrix(rawCounts,colData,formula(~treatment))

# Check DEseq matrix
dds

# Run transformation, blind=T to blind the transformation to the experimental design 
vsd <- varianceStabilizingTransformation(dds, blind=T)

# Pull out transformed data as matrix. rows = genes, columns = individuals
normCounts= assay(vsd)

head(normCounts)

datExpr <- normCounts %>% t() %>% as.data.frame()

# Check samples for excessive missing values or zero variance
gsg = goodSamplesGenes(datExpr, verbose = 3)
# Fails - some genes with zero variance

# Filter out counts that don't pass 
normCounts.goodGenes <- normCounts %>% as.data.frame() %>% mutate(isGood = gsg$goodGenes) %>% filter(isGood) %>% select(-isGood)

datExpr <- normCounts.goodGenes %>% t() %>% as.data.frame()

# Check again
gsg2 = goodSamplesGenes(datExpr, verbose = 3)
# All good


# Cluster samples and look for any obvious outliers
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Generate table of trait data
traitData  <-  sampleInfo %>% 
  mutate(cell_type = str_sub(treatment,1,1),
         serum_type = str_sub(treatment,2,2),
         sex = ifelse(str_detect(sample,'O|C'),0,1), # Female == 0, Male == 1
         hib_cells = ifelse(cell_type == 'H',1,0),
         act_cells = ifelse(cell_type == 'A',1,0),
         hib_serum = ifelse(serum_type == 'H',1,0),
         act_serum = ifelse(serum_type == 'A',1,0),
         glu_serum = ifelse(serum_type == 'G',1,0)) %>% 
  select(-treatment,-cell_type,-serum_type)

head(traitData)
  
# Form a data frame analogous to expression data that will hold the clinical traits.
samples = rownames(datExpr)
traitRows = match(samples, traitData$sample)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage()

# Visualize traits on dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


## Network generation
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType = 'signed hybrid')

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# Construct network w/ power of 8 (*determined above)
cor <- WGCNA::cor # setting cor to be the WGCNA function to avoid errors downstream

# net = blockwiseModules(datExpr, power = 8,
#                        networkType = 'signed hybrid',
#                        TOMType = "signed", minModuleSize = 30, 
#                        mergeCutHeight = 0.25,
#                        numericLabels = T, 
#                        pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "CellCultuerTOM",
#                        verbose = 3)

net = blockwiseModules(datExpr, power = 8, networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE, 
                       #saveTOMFileBase = "fatTOM",
                       maxBlockSize = 14000,
                       verbose = 5)

table(net$colors)

# Plot 
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
   file = "3_analyses/WGCNA-networkConstruction-auto.RData")

# Relating network to external traits

load('3_analyses/WGCNA-networkConstruction-auto.RData')

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Plot
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               plotLegend = F,
               main = paste("Module-trait relationships"))

ME.order <- names(MEs) %>% str_sub(3)

moduleSizes <- as.data.frame(moduleColors) %>% 
  group_by(moduleColors) %>% 
  tally() %>% 
  mutate(moduleColors = factor(moduleColors,levels=rev(ME.order)))


ggplot(moduleSizes,aes(x=n,y=moduleColors,fill=moduleColors)) +
  geom_bar(stat='identity',color='black') +
  geom_text(aes(label=n,x=n+200)) +
  scale_fill_identity() +
  labs(x='Number of Genes',y='') +
  scale_x_continuous(expand=c(0,0),limits = c(0,4000)) +
  theme_classic()

# Get list of genes per module
moduleGeneMembers <- as.data.frame(moduleLabels) %>% 
  mutate(module_color = moduleColors) %>% 
  rownames_to_column('gene')

moduleGeneMembers %>% 
  group_by(module_color) %>% 
  tally()

# write_tsv(moduleGeneMembers,'3_analyses/WGCNA_GeneModuleMembership_05.26.22.tsv')



# Overlap WGCNA modules with DE results -----------------------------------

moduleGeneMembers <- read_tsv('3_analyses/WGCNA_GeneModuleMembership_05.26.22.tsv')

length(unique(moduleGeneMembers$module_color))

sigMods <- c('purple','green','red','black','pink','midnightblue','darkgrey','darkturquoise','brown','darkred','cyan','lightcyan','navajowhite2','blue')

moduleGeneMembers.sig <- moduleGeneMembers %>% filter(module_color %in% sigMods)


de.HA_HH <- read_csv('2_data/de_gene_lists/tt_HA_HH.csv') %>% left_join(moduleGeneMembers.sig,by='gene')
de.AA_HH <- read_csv('2_data/de_gene_lists/tt_AA_HH.csv') %>% left_join(moduleGeneMembers.sig,by='gene')
de.HPD_HH <- read_csv('2_data/de_gene_lists/tt_HPD_HH.csv')%>% left_join(moduleGeneMembers.sig,by='gene')
de.AH_HH <- read_csv('2_data/de_gene_lists/tt_AH_HH.csv')  %>% left_join(moduleGeneMembers.sig,by='gene')
de.AA_AH <- read_csv('2_data/de_gene_lists/tt_AA_AH.csv') %>% left_join(moduleGeneMembers.sig,by='gene')
de.AH_APD <- read_csv('2_data/de_gene_lists/tt_AH_APD.csv')%>% left_join(moduleGeneMembers.sig,by='gene')
de.HA_HPD <- read_csv('2_data/de_gene_lists/tt_HA_HPD.csv')%>% left_join(moduleGeneMembers.sig,by='gene')
de.AA_APD <- read_csv('2_data/de_gene_lists/tt_AA_APD.csv')%>% left_join(moduleGeneMembers.sig,by='gene')

de.all <- de.HA_HH %>% 
  bind_rows(de.AA_HH,de.HPD_HH,de.AH_HH,de.AA_AH,de.AH_APD,de.HA_HPD,de.AA_APD) %>% 
  mutate(comparison = paste(treat1,' vs. ',treat2,sep='')) %>% 
  mutate(comparison = factor(comparison,levels = c('HA vs. HH','AA vs. HH','HG vs. HH','AH vs. HH','AA vs. AH','AH vs. AG','AA vs. AG','HA vs. HG')))

de.all.counts <- de.all %>% 
  group_by(comparison,module_color) %>% 
  tally() 

# ggplot(de.all,aes(x=comparison,fill=module_color)) +
#   geom_bar(stat='count',color='black') +
#   scale_fill_identity() +
#   theme_linedraw()

ggplot(de.all.counts,aes(x=comparison,y=module_color,fill=n)) +
  geom_tile(color='black') + 
  geom_text(data=de.all.counts %>% filter(n<700),aes(label=n),color='white') +
  geom_text(data=de.all.counts %>% filter(n>700),aes(label=n),color='black') +
  scale_fill_continuous_sequential('viridis',rev=F)+
  labs(x='Comparison',y='Module',fill='Number of Genes') +
  theme_linedraw(base_size = 16) +
  scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


##