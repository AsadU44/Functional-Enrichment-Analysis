library(clusterProfiler)
library(DOSE)
library(R.utils)
library(gtools)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)
library(tibble)
library(ggupset)
library(DOSE)
library(AnnotationDbi)
library(pathview)

setwd('E:/R-Programming-Practices/Functional Enrichment/Gene Ontology Terms')
### Read data
GENEs<- data.frame(read.csv('GBM.csv'))

### Convert Symbol to ENTREZID
Converted <- bitr(unique(GENEs$Genes), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db) 
head(Converted)


### Prepare data for Functional Enrichment Analysis
Gene_List <- Converted$ENTREZID
Gene_List = sort(Gene_List, decreasing = TRUE) 
head(Gene_List)
#N.B. Sorting to deacreasing order is required for enrichment analysis

### Overrepresentation (ORA) Analysis altogether
All_result <- enrichGO(gene     = Gene_List,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.01)

#### Visualization
barplot(All_result, split = "ONTOLOGY")+facet_grid(ONTOLOGY~., scale = "free")+
                                                   ggtitle("Barplot for ORA")
### Biological processes gene ontlogy terms
BP_result <- enrichGO(gene     = Gene_List,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP", #Should be changed during MF and CC
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.01)
#### Visualization
dotplot(BP_result, showCategory=15) + ggtitle("Dotplot for BP Go Terms")

### Integrate log2FC for ridge plot
#In this step we need prepare our data similar to geneList (package='DOSE')
#Let's read new data and omit any NA values
DEGs<-read.csv('DEGs.csv')
DEGs<- na.omit(DEGs)

#### Convert gene symbols to ENTREZID
Converted <- bitr(unique(DEGs$Genes), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db) 
head(Converted)

#### Prepare data
gene_list<- DEGs$logFC
names(gene_list)<- Converted$ENTREZID
gene_list<- sort(gene_list, decreasing = T)

#### Let's perform disease ontology (DO) analysis
DO<- gseDO(gene_list, pvalueCutoff = 1)
ridgeplot(DO, showCategory = 10)

### Disease gene network (DGN) Analysis
#### Prepare data
dgn_data<-names(gene_list)[abs(gene_list) > 2]
dgn_data<-sort(dgn_data, decreasing = T)

#### Perform DGN Analysis
dgn_res <- enrichDGN(dgn_data, pvalueCutoff = 0.05) 
head(dgn_res)

#### Convert gene id to name in the enrichment object
dgn_res <- setReadable(dgn_res, 'org.Hs.eg.db', 'ENTREZID')

#### Plot
Plot1<- heatplot(dgn_res, showCategory = 15) + ggtitle("Heatmap of DGN Analysis")
Plot1
Plot2<- heatplot(dgn_res, showCategory = 15, foldChange=gene_list) + 
        ggtitle("Heatmap of DGN Analysis") 
Plot2

#### Combine two plots together
cowplot::plot_grid(Plot1, Plot2, ncol=1, labels=LETTERS[1:2])

#### Tree plot
tree_data<- pairwise_termsim(dgn_res)
Tree_plot<-treeplot(tree_data)
Tree_plot

#### Cnet plot
Cnet_plot <- cnetplot(dgn_res, foldChange=gene_list, circular = TRUE,
                      showCategory=10, colorEdge = TRUE)
Cnet_plot

#N.B. Cnetplot, Tree plot and Heatplot work best for small number of genes. 

#Reference: 'http://yulab-smu.top/clusterProfiler-book/index.html'