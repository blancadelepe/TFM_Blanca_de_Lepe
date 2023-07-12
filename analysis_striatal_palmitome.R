#################################################################
##   Analysis of the striatal palmitome from zdhcc15-KO mice   ##
##                                                             ##
##                Author: Blanca de Lepe López                 ##
#################################################################

# Load the data and rearrange it
intensity.raw.data <- read.table(file = "data/intensidad_proteinas.tsv",
                                 header = T, sep="\t", as.is = T)

prot.ids <- intensity.raw.data$ProtID
rownames(intensity.raw.data) <- prot.ids
intensity.raw.data = intensity.raw.data[,-1]
intensity.raw.data[intensity.raw.data==0] = NA

################################################################################
################################################################################


# Qualitative palmitoylated proteins are extracted, which will be those that 
# have all the values in one of the groups, and none of the values in the other 
# group.
# A function is created to extract the proteins that have no value in the KO 
# group and all of them in the WT group.

detectar_filas_knockout = function(df) {
  filas_con_nas = apply(df[, 1:3], 1, function(x) sum(is.na(x)) == 3) 
  filas_sin_na = apply(df[, (ncol(df)-2):ncol(df)], 1, function(x) sum(is.na(x)) == 0) 
  filas_validas = filas_con_nas & filas_sin_na 
  return(df[filas_validas, ])
}

ko_qualitative = detectar_filas_knockout(intensity.raw.data)
ko_qualitative_vector = row.names(ko_qualitative)

# A function is created to extract the proteins that have no value in the WT 
# group and all of them in the KO group.

detectar_filas_wild = function(df) {
  filas_con_nas = apply(df[, (ncol(df)-2):ncol(df)], 1, function(x) sum(is.na(x)) == 3)
  filas_sin_na = apply(df[, 1:3], 1, function(x) sum(is.na(x)) == 0)
  filas_validas = filas_con_nas & filas_sin_na
  return(df[filas_validas, ])
}


wt_qualitative = detectar_filas_wild(intensity.raw.data)
wt_qualitative_vector = row.names(wt_qualitative)

qualitative.proteins = c(ko_qualitative_vector, wt_qualitative_vector)
write.table(x = qualitative.proteins,file = "tables/qualitative_proteins.txt",quote = F,
            row.names = F, sep = "\t", )
qualitative.proteins.gene = read.table(file = "data/qualitative.proteins.to.genes.txt",
                                       header = T, sep = "\t")
qualitative.proteins.gene = qualitative.proteins.gene[,2]

# Value imputation (replacement of NAs by the mean of the rows)

intensity.raw.data <- lapply(intensity.raw.data, as.numeric)
intensity.raw.data <- as.data.frame(intensity.raw.data)
rownames(intensity.raw.data) <- prot.ids

# Divide the dataframe in 2 in order to be able to study each sample 
# independently 
KO = intensity.raw.data[,1:3]; WT = intensity.raw.data[,4:6]

# The rows of each group in which more than 1 NA appears are deleted.
KO = KO[rowSums(is.na(KO)) < 2, ]
WT = WT[rowSums(is.na(WT)) < 2, ]

KO$mean.expr <- rowMeans(KO, na.rm = T)
for (i in 1:nrow(KO)){
  for (j in 1:ncol(KO)){
    if (is.na(KO[i,j])==TRUE){
      KO[i,j] <- KO[i,4]
    }
  }
}

WT$mean.expr <- rowMeans(WT, na.rm = T)
for (i in 1:nrow(WT)){
  for (j in 1:ncol(WT)){
    if (is.na(WT[i,j])==TRUE){
      WT[i,j] <- WT[i,4]
    }
  }
}


interseccion = intersect(rownames(KO[,1:3]),rownames(WT[,1:3]))

intensity.raw.data = intensity.raw.data[interseccion,]
dim(intensity.raw.data)
prot.ids.modified <- row.names(intensity.raw.data)

write.table(x = intensity.raw.data,file = "tables/intensity.raw.data.modified.tsv",quote = F,
            row.names = F, sep = "\t")

################################################################################
################################################################################

library(NormalyzerDE)
library(MetBrewer)
library(tidyverse)


# Visualisation of the original data before normalization

boxplot(intensity.raw.data, outline = F, main="Boxplot preprocessed data",cex.main=1.2, 
        col=met.brewer(n=6, name="Hiroshige"), las=2)
title(ylab = "LFQ intensity (not normalized)", xlab="Samples", line = 4)

intensity.raw.data %>% 
  gather(key="Samples", value="Val") %>%
  ggplot( aes(x=Samples, y=Val, fill=Samples)) +
  geom_violin(width=0.9) +
  geom_boxplot(width=0.1, color="black", alpha=0.2) + labs(x='Samples',
                                                           y='LFQ intensity (not normalized)') + theme_classic()

# Data are prepared for normalization 

design <- data.frame(sample=colnames(intensity.raw.data),
                     group=c(rep("KO",3),rep("WT",3)))
design

write.table(x = design,file = "tables/intensity.data.design.tsv",quote = F,row.names = F,
            sep = "\t")

# Normalization with the normalyzer function

normalyzer(jobName = "normalization.intensity", designPath = "tables/intensity.data.design.tsv",
           dataPath = "tables/intensity.raw.data.modified.tsv" ,outputDir = ".")

# Reading normalised data

intensity.normalized.data <- read.table(file="normalization.intensity//Quantile-normalized.txt",
                                    header = T,as.is=T,sep="\t")

row.names(intensity.normalized.data) <- prot.ids.modified
head(intensity.normalized.data)


# Visualisation of the normalized data 

boxplot(intensity.normalized.data, outline = F, main="Boxplot normalized data",cex.main=1.2, 
        col=met.brewer(n=6, name="Hiroshige"), las=2, cex.lab= 1, 
        ylab = "LFQ intensity (normalized)")


intensity.normalized.data %>% 
 gather(key="Samples", value="Val") %>%
 ggplot( aes(x=Samples, y=Val, fill=Samples)) +
 geom_violin(width=0.9) +
 geom_boxplot(width=0.1, color="black", alpha=0.2) + labs(x='Samples',
                                                           y='LFQ intensity (normalized)') + theme_classic()

################################################################################
################################################################################

# Here we perform a PCA and a dendrogram

library(FactoMineR)
library(factoextra)

pca.normalized.data <- data.frame(t(intensity.normalized.data))

res.pca <- PCA(pca.normalized.data, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 2)  
fviz_dend(res.hcpc,k=2,
          cex = 0.75,                       
          palette = "jco",              
          rect = TRUE, rect_fill = TRUE, 
          rect_border = "jco",           
          type="rectangle",
          labels_track_height = 1400)

fviz_pca_ind(res.pca, col.ind = c("KO1","KO2","KO3","WT1","WT2","WT3"), 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE)

################################################################################
################################################################################

## Differentially palmitoylated protein analysis with limma

library(limma)

# The experimental design is specified.

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2)))
colnames(experimental.design) <- c("KO","WT")


# The estimate of the intensity levels of each protein is fitted to a linear 
# model taking into account the experimental design.

linear.fit <- lmFit(intensity.normalized.data, experimental.design)

# Contrast is specified.

contrast.matrix <- makeContrasts(KO-WT,levels=c("KO","WT"))

# The fold-change and corresponding p-values are calculated for each protein in
# the specified contrast using the *constrasts.fit* and *eBayes* functions.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

# The topTable function performs a difference analysis of the palmitoylation 
# levels of the proteins, and returns a table with information about the 
# proteins and their intensity in the different groups.

WT.KO <- topTable(contrast.results, number=1236,coef=1,sort.by="logFC")

log.fold.change <- WT.KO$logFC
q.value <- WT.KO$adj.P.Val
p.value <- WT.KO$P.Value
prot.ids.WT.KO <- rownames(WT.KO)
names(log.fold.change) <- prot.ids.WT.KO
names(q.value) <- prot.ids.WT.KO
names(p.value) <- prot.ids.WT.KO


# Protein filtering using logFC of 2, 1.6, 1 and 0.6 and a q-value of 0.05.

q.more.palmitoylated.2 <- prot.ids.WT.KO[log.fold.change > 2 & q.value < 0.05]
q.less.palmitoylated.2 <- prot.ids.WT.KO[log.fold.change < -2 & q.value < 0.05]

q.more.palmitoylated.1.6 <- prot.ids.WT.KO[log.fold.change > 1.6 & q.value < 0.05]
q.less.palmitoylated.1.6 <- prot.ids.WT.KO[log.fold.change < -1.6 & q.value < 0.05]

q.more.palmitoylated.1 <- prot.ids.WT.KO[log.fold.change > 1 & q.value < 0.05]
q.less.palmitoylated.1 <- prot.ids.WT.KO[log.fold.change < -1 & q.value < 0.05]

q.more.palmitoylated.0.6 <- prot.ids.WT.KO[log.fold.change > 0.6 & q.value < 0.05]
q.less.palmitoylated.0.6 <- prot.ids.WT.KO[log.fold.change < -0.6 & q.value < 0.05]

length(q.more.palmitoylated.2)
length(q.less.palmitoylated.2)

length(q.more.palmitoylated.1.6)
length(q.less.palmitoylated.1.6)

length(q.more.palmitoylated.1)
length(q.less.palmitoylated.1)

length(q.more.palmitoylated.0.6)
length(q.less.palmitoylated.0.6)

# Only one protein was found to have a significantly lower intensity with a 
# logFC of 0.6 and a q.value of 0.5.



# Protein filtering using logFC of 2, 1.6, 1 and 0.6 and a p-value of 0.05.


more.palmitoylated.2 <- prot.ids.WT.KO[log.fold.change > 2 & p.value < 0.05]
less.palmitoylated.2 <- prot.ids.WT.KO[log.fold.change < -2 & p.value < 0.05]

more.palmitoylated.1.6 <- prot.ids.WT.KO[log.fold.change > 1.6 & p.value < 0.05]
less.palmitoylated.1.6 <- prot.ids.WT.KO[log.fold.change < -1.6 & p.value < 0.05]

more.palmitoylated.1 <- prot.ids.WT.KO[log.fold.change > 1 & p.value < 0.05]
less.palmitoylated.1 <- prot.ids.WT.KO[log.fold.change < -1 & p.value < 0.05]

more.palmitoylated.0.6 <- prot.ids.WT.KO[log.fold.change > 0.6 & p.value < 0.05]
less.palmitoylated.0.6 <- prot.ids.WT.KO[log.fold.change < -0.6 & p.value < 0.05]



length(more.palmitoylated.2)
length(less.palmitoylated.2)

length(more.palmitoylated.1.6)
length(less.palmitoylated.1.6)

length(more.palmitoylated.1)
length(less.palmitoylated.1)

length(more.palmitoylated.0.6)
length(less.palmitoylated.0.6)


write.table(more.palmitoylated.0.6, file="tables/more.palmitoylated.0.6.tsv", 
            sep="\t", quote = FALSE, col.names = T, row.names = FALSE)
write.table(less.palmitoylated.0.6, file="tables/less.palmitoylated.0.6.tsv", 
            sep="\t", quote = FALSE, col.names = T, row.names = FALSE)

# From now on, we will work with the 34 proteins with higher levels of 
# palmitoylation and the 40 proteins with lower levels of palmitoylation that 
# have been obtained by filtering by a logFC of 0.6 and a p-value of 0.05.

################################################################################
################################################################################

# We extract the list of differentially palmitoylated proteins with their 
# respective logFC and p-values.

library(dplyr)
library(tibble)


WT.KO = WT.KO %>%
  rownames_to_column(var = "Protein")

table.more = WT.KO %>%
  filter(Protein %in% more.palmitoylated.0.6) %>%
  select(logFC, P.Value)

row.names(table.more) = more.palmitoylated.0.6

write.table(table.more, file="tables/more.logfc.pvalue.tsv", 
            sep="\t", quote = FALSE, col.names = T, row.names = FALSE)


table.less = WT.KO %>%
  filter(Protein %in% less.palmitoylated.0.6) %>%
  select(logFC, P.Value)

row.names(table.less) = less.palmitoylated.0.6

write.table(table.less, file="tables/less.logfc.pvalue.tsv", 
            sep="\t", quote = FALSE, col.names = T, row.names = FALSE)

################################################################################
################################################################################

# Volcano plots for visualization of differentially palmitoylated proteins

log.p.val <- -log10(p.value)

plot(log.fold.change,log.p.val,pch=19,col="grey",cex=0.8,
     xlim=c(-12,12),ylim = c(0,4), 
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1)

points(x = log.fold.change[more.palmitoylated.0.6],
       y = log.p.val[more.palmitoylated.0.6],col="red",cex=0.8,pch=19)
points(x = log.fold.change[less.palmitoylated.0.6],
       y = log.p.val[less.palmitoylated.0.6],col="blue",cex=0.8,pch=19)
title(main = "Differential palmitoylated proteins")
legend("topleft", legend = c("More palmitoylated", "Less palmitoylated"), 
       fill = c("red", "blue"), cex = 0.5)

################################################################################
################################################################################

## Functional enrichment analysis of Gene Ontology terms

library(clusterProfiler)
library(org.Mm.eg.db)

# More palmitoylated proteins functional enrichment by ontology

more.palmitoylated.enrich.go.cc <- enrichGO(gene = more.palmitoylated.0.6,
                                     OrgDb         = "org.Mm.eg.db",
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable      = T,
                                     keyType = "UNIPROT")

barplot(more.palmitoylated.enrich.go.cc, showCategory = 8)

more.palmitoylated.enrich.go.mf <- enrichGO(gene = more.palmitoylated.0.6,
                                         OrgDb         = "org.Mm.eg.db",
                                         ont           = "MF",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = T,
                                         keyType = "UNIPROT")

# No enriched terms found for molecular function


more.palmitoylated.enrich.go.bp <- enrichGO(gene = more.palmitoylated.0.6,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "UNIPROT")

# No enriched terms found for biological process

# Less palmitoylated proteins functional enrichment by ontology

less.palmitoylated.enrich.go.cc <- enrichGO(gene = less.palmitoylated.0.6,
                                         OrgDb         = "org.Mm.eg.db",
                                         ont           = "CC",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = T,
                                         keyType = "UNIPROT")

barplot(less.palmitoylated.enrich.go.cc,showCategory = 10)


less.palmitoylated.enrich.go.mf <- enrichGO(gene = less.palmitoylated.0.6,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "MF",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "UNIPROT")

barplot(less.palmitoylated.enrich.go.mf,showCategory = 10)


less.palmitoylated.enrich.go.bp <- enrichGO(gene = less.palmitoylated.0.6,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "UNIPROT")

barplot(less.palmitoylated.enrich.go.bp,showCategory = 10)


# Less and more palmitoylated proteins functional enrichment with all ontologies 

less.palmitoylated.enrich.go.all <- enrichGO(gene = less.palmitoylated.0.6,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "ALL",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "UNIPROT")

barplot(less.palmitoylated.enrich.go.all,showCategory = 10)


more.palm.gene.symbol <- read.table(file = "data/more.palm.to.genes.txt", header = T, sep = "\t")
more.palm.gene.symbol <- more.palm.gene.symbol[,2]

less.palm.gene.symbol <- read.table(file = "data/less.palm.to.genes.txt", header = T, sep = "\t")
less.palm.gene.symbol <- less.palm.gene.symbol[,2]

more.palmitoylated.enrich.go.all <- enrichGO(gene = more.palm.gene.symbol,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "ALL",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "SYMBOL")

barplot(more.palmitoylated.enrich.go.all, showCategory = 15)
cnetplot(more.palmitoylated.enrich.go.all,showCategory = 15, cex.params = list(category_label = 0.7), )
heatplot(more.palmitoylated.enrich.go.all, showCategory = 15)


less.palmitoylated.enrich.go.all <- enrichGO(gene = less.palm.gene.symbol,
                                            OrgDb         = "org.Mm.eg.db",
                                            ont           = "ALL",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.05,
                                            readable      = T,
                                            keyType = "SYMBOL")

barplot(less.palmitoylated.enrich.go.all, showCategory = 15)
cnetplot(less.palmitoylated.enrich.go.all,showCategory = 15, cex.params = list(category_label = 0.7))
heatplot(less.palmitoylated.enrich.go.all, showCategory = 15)


# Differentially palmitoylated proteins functional enrichment

diff.palm.gene.symbol <- read.table(file = "data/diff.proteins.to.genes.txt", header = T, sep = "\t")
diff.palm.gene.symbol <- diff.palm.gene.symbol[,2]

diff.palmitoylated.enrich.go.all <- enrichGO(gene = diff.palm.gene.symbol,
                                             OrgDb         = "org.Mm.eg.db",
                                             ont           = "ALL",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.05,
                                             readable      = T,
                                             keyType = "SYMBOL")

barplot(diff.palmitoylated.enrich.go.all, showCategory = 15)
cnetplot(diff.palmitoylated.enrich.go.all,showCategory = 12, cex.params = list(category_label = 0.7))
heatplot(diff.palmitoylated.enrich.go.all, showCategory = 20)


# Qualitative palmitoylated proteins functional enrichment

qualitative.proteins.enrich.go <- enrichGO(gene = qualitative.proteins.gene,
                                             OrgDb         = "org.Mm.eg.db",
                                             ont           = "ALL",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.05,
                                             readable      = T,
                                             keyType = "SYMBOL")

barplot(qualitative.proteins.enrich.go, showCategory = 20)

