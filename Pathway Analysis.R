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

#Ser wd and read data
setwd('E:/R-Programming-Practices/Functional Enrichment/Pathway Analysis')
data= data.frame(read.csv('GBM.csv'))

#Convert gene symbol to ENTREZID
Converted <- bitr(unique(data$Genes), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db) 
head(Converted)

#Prepare gene list for KEGG analysis
gene_list <- Converted$ENTREZID
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)

#Perform KEGG analysis
kegg <- enrichKEGG(gene         = gene_list,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kegg)

#Visualization
#Dotplot
dotplot(kegg, showCategory=15) + ggtitle("Barplot for KEGG Analysis")

#Barplot
barplot(kegg, showCategory=15) + ggtitle("Barplot for KEGG Analysis")

#Check KEGG enrichment of the input gene list online in one pathway
hsa04024<- pathview(gene.data  = gene_list,
                    pathway.id = "hsa04024",
                    species    = "hsa")


#Check KEGG enrichment of the input gene list offline (wd) in one pathway
hsa04110 <- pathview(gene.data  = gene_list,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     )
                     
#Integrate log2FC in KEGG pathway
DEGs=read.csv('DEGs.csv')

Converted <- bitr(unique(DEGs$Genes), fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db) 
head(Converted)

#Prepare gene list for KEGG analysis
Gene_List <- Converted$ENTREZID
Gene_List = sort(Gene_List, decreasing = TRUE)
head(Gene_List)

#Create object
logFC<- DEGs[,2]
names(logFC)<-Gene_List
mypathway<-"hsa04110"

#View pathway
pathview(gene.data=logFC,species="hsa",
         pathway=mypathway, limit = list(gene=3, cpd=1))

#Wiki_Pathway Analysis
WP=enrichWP(gene_list, organism = 'Homo sapiens')
ridgeplot(WP, showCategory = 15) + ggtitle("Barplot for WikiPathway Analysis")
