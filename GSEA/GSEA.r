library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(DOSE)

detach(package:clusterProfiler, unload = T)
data<- read.csv('DEGs.csv')

data<- data[order(-data$log2FoldChange),]

gene_list<- data$log2FoldChange
names(gene_list)<- data$X

gse<- gseGO(gene_list,
            keyType = 'ENSEMBL',
            OrgDb = 'org.Hs.eg.db',
            ont = 'BP',
            pvalueCutoff = 1
            )

as.data.frame(gse)

#General plot
plot<- gseaplot2(gse, geneSetID = 2, title = gse$Description[2]) 
plot

#Multiple entichment item
plot2<- gseaplot2(gse, geneSetID = 1:5, color = c('#D43F3A', '#EEA236', '#5CB85C',
                                                  '#46B8DA', '#9632B8')) 
plot2

#Add p value table
plot3<- gseaplot2(gse, geneSetID = 1:3, color = c('#D43F3A', '#EEA236', '#5CB85C'
                 ), pvalue_table = T) 
plot3

#Upset plot
upsetplot(gse)

#For Diseae Ontology Term
Converted <- bitr(unique(data$X), fromType = "ENSEMBL", toType = c(
                    "ENTREZID"), OrgDb = org.Hs.eg.db)
Gene_List<- Converted$ENTREZID
Gene_List = sort(Gene_List, decreasing =TRUE) 
head(Gene_List)

DGN<-enrichDO(Gene_List)

Plot4<- upsetplot(DGN)
Plot4