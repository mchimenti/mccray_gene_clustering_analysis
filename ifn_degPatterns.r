##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")
#biocLite("COMBINE-lab/wasabi")       #install wasabi the first time

#negative binomial GLM
library(DESeq2)

#annotation
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library(tidyverse)
#library(pcaExplorer)
#pathway
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/mccray/project_cf_rsv_sinnbrajesh_human_nov2017/")

#######################################
## FeatureCounts > DESeq2 > PCAExplorer
#######################################

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## setup paths and create count matrix and metadata table; create DESeq object

matrix_path <- "/Users/mchimenti/iihg/RNA_seq/mccray/project_cf_rsv_sinnbrajesh_human_nov2017/2017-11-10_bcbio_run_dir/annotated_combined.counts"
count_mat <- read.table(matrix_path, header=TRUE, row.names = 'id')

annot <- count_mat[49]
count_mat <- count_mat[,1:48]
colnames(count_mat)[1] <- "RSV_NCF2_0hrs"  ## for this exp. MeV == RSV.  
count_mat <- count_mat[,c(2:31,1,32:48)]  #put the columns in order for ease of metadata addition
count_mat <- as.matrix(count_mat)

CF <- paste(rep("CF", 24))
NCF <- paste(rep("NCF", 24))
condition <- c(CF, NCF)

times <- c("h0","h920","h12","h24","h48","h72")  # 'h920' stands in for 'h120' to get R to sort the factor correctly
timecourse <- rep(times, 8)

colData <- as.data.frame(condition, row.names = colnames(count_mat))
colData['time'] <- timecourse
all(rownames(colData) %in% colnames(count_mat))

colData$group <- factor(paste0(colData$condition, colData$time))


#########
##DE Analysis combining two factors into "group" factor for 2-way comparo

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = colData,
                              design = ~ group)

## prefiltering and re-leveling

dds <- dds[ rowSums(counts(dds)) > 5, ]
dds$group <- relevel(dds$group, ref='NCFh0')

## Differential expression analysis with contrasts

dds <- DESeq(dds)

anno <- get_annotation(dds, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rld <- rlog(dds, blind=FALSE)


########
## Clustering on IFN gene expression only
########

## get IFN gene set
ifn_genes <- read.csv('~/iihg/RNA_seq/mccray/project_gene_clustering_cf_rsv_sept2018/typeI_interferon_inducible_genes.csv', header=FALSE)
names(ifn_genes) <- "gene_name"
ifn_genes <- left_join(ifn_genes, anno, on = 'gene_name')
ifn_genes <- na.omit(ifn_genes) # drop not found gene names 
ifn_genes <- ifn_genes[ifn_genes$gene_id %in% row.names(counts(dds)), ]  # drop transcripts not in the SO 

## clustering 
library(DEGreport)

## clustering on whole experiment (CF and NCF)
ifn_mat <- assay(rld)[ifn_genes$gene_id,]
clusters <- degPatterns(ifn_mat, metadata = as.data.frame(colData(dds)), time = 'time', col = 'condition',
                        reduce = TRUE, cutoff = 0.85, minc = 20)


## clustering on just CF or NCF 
cluster_cf <- degPatterns(ifn_mat[,1:24], metadata = as.data.frame(colData(dds)[1:24,]), time = 'time')
cluster_ncf <- degPatterns(ifn_mat[,25:48], metadata = as.data.frame(colData(dds)[25:48,]), time = 'time')

###########
## Clustering on all DE genes btw CF and NCF at 24 hours
###########

res_CF_24 <- results(dds, contrast = c("group", "CFh24", "NCFh24"))
res_CF_24 <- na.omit(res_CF_24)  #drop NA rows
res_CF_24_sig <- res_CF_24[res_CF_24$padj < 0.001 & res_CF_24$baseMean > 5.0,]
res_CF_24_ord <- res_CF_24_sig[order(res_CF_24_sig$padj),]


CF_24_mat <- assay(rld)[row.names(res_CF_24_ord), ]
clusters_CF_24 <- degPatterns(CF_24_mat, metadata = as.data.frame(colData(dds)), time = 'time', col = 'condition', minc = 15)

res_CF_72 <- results(dds, contrast = c("group", "CFh72", "NCFh72"))
res_CF_72 <- na.omit(res_CF_72)  #drop NA rows
res_CF_72_sig <- res_CF_72[res_CF_72$padj < 0.05 & res_CF_72$baseMean > 5.0,]
res_CF_72_ord <- res_CF_72_sig[order(res_CF_72_sig$padj),]


CF_72_mat <- assay(rld)[row.names(res_CF_72_ord), ]
clusters_CF_72 <- degPatterns(CF_72_mat, metadata = as.data.frame(colData(dds)), time = 'time', col = 'condition', minc = 15)


############
## Cluster profiling 
############

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")

library(clusterProfiler)

############
## 24hrs DE genes
############
clust24 <- unique(clusters_CF_24$df$cluster)
for(i in clust24){
  group <- filter(clusters_CF_24$df, cluster %in% as.character(i))
  group.df <- bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = group.df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  readable = TRUE)
  print("Enriched GO terms for cluster")
  print(i)
  print(head(ego))
}

for(i in clust24){
  group <- filter(clusters_CF_24$df, cluster %in% as.character(i))
  group.df <- suppressMessages(bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  ego <- enrichKEGG(gene = group.df$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
  print(paste0("Enriched KEGG terms for cluster #:", i))
  print(head(ego))
}


###############
## Clustering at DE genes at 72 hours
###############

clust <- unique(clusters_CF_72$df$cluster)
for(i in clust){
  group <- filter(clusters_CF_72$df, cluster %in% as.character(i))
  group.df <- bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = group.df$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  readable = TRUE)
  print("Enriched GO terms for cluster")
  print(i)
  print(head(ego))
}

for(i in clust){
  group <- filter(clusters_CF_72$df, cluster %in% as.character(i))
  group.df <- suppressMessages(bitr(group$genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  ego <- enrichKEGG(gene = group.df$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)
  print(paste0("Enriched KEGG terms for cluster #:", i))
  print(head(ego))
}


