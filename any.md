``` r
mat <- read.csv("rna_tpm_log2.csv", row.names = 1)
```

Generating data frame for visualization

``` r
test <- mat["PDK1",]
test <- stack(as.data.frame(test))
test$ind <- colnames(mat)
test$type <- ifelse(test$ind %in% c("MKN1", "Hs746T", "SNU484", "SNU668", "YCC16", "YCC11", "SK4", "SNU1750"), "SEM", "nSEM")
colnames(test) <- c("Expression_Level", "Cell", "Type")
write.csv(test, "PDK1.csv")


library(ggplot2)
library(ggpubr)
```

Compare between groups (Type) + wilcox-rank sum test
``` r
p <- ggboxplot(test, x = "Type", y = "Expression_Level", ylab = "Expression_Level",fill = "Type", 
               title = "HMGCR") + 
  theme(axis.text.y=element_text(size=30),
        axis.text.x = element_text(angle = , vjust = 0.5, hjust=1, size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        plot.title = element_text(size = 40, face = "bold"))
a <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               plot.background = element_rect(fill = "transparent", color = "white"),
               panel.background = element_rect(fill = "transparent", color = "white"))

a + stat_compare_means( aes(label = ..p.signif..),) + stat_compare_means(size = 5)
```


##################Visualization of all elements - In this case, Descending Order######################
``` r
p <- ggplot(test) + 
  geom_boxplot(aes(x= reorder(Cell, -Expression_Level), y = Expression_Level), ) +
  ggtitle("HMGCR") +
  theme(axis.text.y=element_text(size=30),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40, face = "bold"),
        plot.title = element_text(size = 40, face = "bold"))

a <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               plot.background = element_rect(fill = "transparent", color = "white"),
               panel.background = element_rect(fill = "transparent", color = "white"))
```

Heatmap between two groups
``` r
coldata <- data.frame(matrix(0, ncol = 2, nrow = ncol(mat)))
colnames(coldata) <- c("Cell", "Type")
coldata$Cell <- colnames(mat)
coldata$Type <- ifelse(coldata$Cell %in% c("MKN1", "Hs746T", "SNU484", "SNU668", "YCC16", "YCC11", "SK4", "SNU1750"), "SEM", "nSEM")

colours <- list('Type' = c('SEM' = 'red2', 'nSEM' = 'royalblue'))
colAnn <- HeatmapAnnotation(df = coldata$Type,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),
                            simple_anno_size = unit(0.1, "cm"))
split <- factor(paste0("", coldata$Type))
library(circlize)
library(ComplexHeatmap)

mat_scaled <- t(scale(t(mat)))
#cairo_pdf("/Users/ryanmachine/Datadriven/kmeans(5).pdf", height = 6, width = 8)
Heatmap(as.matrix(mat_scaled), name = "Z-score", 
        top_annotation = colAnn, column_split = coldata$Type,
        row_names_gp = gpar(fontsize = 11), column_names_gp = gpar(fontsize = 11), show_row_names = TRUE, 
        show_column_dend = FALSE, use_raster = TRUE, raster_by_magick = TRUE, raster_quality = 2)
```


Modifying Drug target dataframe
-> row : drugname, column: target genes

1. If targets in the dataframe are not in common with provided target, mark them with NA
``` r
AGS <- target
for (j in 1:nrow(AGS)){
  for (i in 3:ncol(AGS)){
    if(!(AGS[j, i] %in% protein$YCC7)){
      AGS[j, i] <- NA
    }
  }
}
```
2. Mark all empty space with NA
``` r
AGS1 <- as.data.frame(apply(AGS, 2, function(x) gsub("^$|^ $", NA, x)))
```

3. All of element to the left (targets left, NAs right)
``` r
new <- as.data.frame(t(apply(AGS1,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
```

4. Get rid of columns with all NAs
``` r
new <- new[!sapply(new, function(x) all(is.na(x)))]
```

5. Convert all NAs to empty space
``` r
new[is.na(new)] <- ""
```

6. Count gene numbers in the dataframe
``` r
newfunc <- new[,-c(1:2)]
solution<-as.data.frame(table(unlist(newfunc)))
solution <- solution[-1,]
```

Converting and pre-processing Microarra raw data to analysis-ready matrix
``` r
cell <- read.delim("gst_i497_cellline/GST_cellLine/TableControl_108_process-3_norm.txt", row.names = 1)
id <- read.delim("ID.txt", row.names = 1)
newcell <- cell[row.names(id), ]
newcell <- newcell[!is.na(rowSums(newcell)),]
id <- id[row.names(newcell),]
id <- subset(id, select = ILMN_Gene)

library(data.table)
newcell$Gene <- id$ILMN_Gene
result <- as.data.frame(do.call(rbind,lapply(lapply(split(newcell,newcell$Gene),`[`,1:ncol(newcell)-1),colMeans)))
converted <- log2(result + 1)
```

K-means clustering and Silhouette method
``` r
library(cluster)
avg_sil <- function(k, data) {
  km.res <- kmeans(data, centers = k)
  ss <- silhouette(km.res$cluster, dist(data))
  avgSil <- mean(ss[, 3])
  return(avgSil)
}
kClusters <- 2:40
resultForEachK <- data.frame(k = kClusters, silAvg = rep(NA, length(kClusters)))
for(i in 1:length(kClusters)){
       resultForEachK$silAvg[i] <- avg_sil(kClusters[i], dataconvert)
}

plot(resultForEachK$k, resultForEachK$silAvg, type = "b", pch = 19, frame = FALSE, xlab = "Number of clusters K",ylab = "Average Silhouettes")

gene_cluster <- kmeans(dataconvert, 30)
genedata <- data.frame(matrix(0, ncol = 2, nrow = nrow(dataconvert)))
genedata$X1 <- rownames(dataconvert)
genedata$X2 <- gene_cluster$cluster
```


NMF-clustering
``` r
library(NMF)
library(nnTensor)
library(ConsensusClusterPlus)

cc.res <- ConsensusClusterPlus(as.matrix(datayonsei), maxK = 5, plot = "png", title = "ALL")
ALL$cc.group <- cc.res[[2]]$consensusClass

nmf.res <- nmf(datayonsei, rank = 2:7)

#NMF clustering
nout <- nmf(datayonsei, 4)
#matrices <- nout@fit@W

pdf("/Users/ryanmachine/factor.pdf")
consensusmap(nmf.res$fit)
basismap(nmf.res$fit)
plot(nmf.res)
dev.off()

breaks = seq(-2, 2, 0.01)
nmf.group <- as.data.frame(predict(nmf.res$fit[[3]]))
colnames(nmf.group) <- "Cluster"
```
TCGA GC RNA-seq data download
``` r
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(sesame)
library(sesameData)
library(EDA)
library(EDASeq)
library(Cairo)
library(biomaRt)
sesameDataCacheAll()

#All data downloaded are Hg38 based harmonized data (Not Legacy)

#Gene expression data
query.met.stad <- GDCquery(
  project = "TCGA-STAD", 
  legacy = FALSE,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)
GDCdownload(query.met.stad, method = "api")

stad.exp <- TCGAbiolinks::GDCprepare(
  query = query.met.stad
)

stadmatrix <- assay(stad.exp)
stad1 <- colData(stad.exp)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = stad.exp, geneInfo = geneInfoHT, method = "gcContent")
write.csv(dataNorm, "TCGA_RNAseq-Normalized.csv")

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

#Download TCGA STAD Clinical data
stad.clin <- GDCquery_clinic(project = "TCGA-STAD", type = "Clinical")
#Distribute clinical data into tumor vs matched normal
clin.tumor <- subset(stad.clin, submitter_id %in% sample_info_tumor$patient)
ciin.normal <- subset(stad.clin, submitter_id %in% sample_info_solid_tissue_normal$patient)

#Extract sample information
sample_info <- as.data.frame(colData(stad.exp))
save(sample_info, file = "TCGA_STAD_sampleinfo.rda")
load("TCGA_STAD_sampleinfo.rda")

colnames(dataFilt) <- substr(colnames(dataFilt), 1, 16)
stad.clin_matched <- subset(stad.clin, submitter_id %in% sample_info$patient)


write.csv(samplesNT, "TCGA-STAD-normal.csv")
write.csv(samplesTP, "TCGA-STAD-tumor.csv")
```
