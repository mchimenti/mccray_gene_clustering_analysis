###  Use TcGSA method to cluster RSV infection timecourse data to find gene sets that 
###  change between CF and nCF types 

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
library(gage)
library(ggplot2)
#timecourse clustering
library(TcGSA)
library(clusterProfiler)

setwd("~/iihg/RNA_seq/mccray/project_gene_clustering_cf_rsv_sept2018/")

#######################################
## FeatureCounts > DESeq2 > PCAExplorer
#######################################

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

## Create the DESEQ object and then the RLD object 
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = colData,
                              design = ~ group)

dds <- dds[ rowSums(counts(dds)) > 5, ]
dds$group <- relevel(dds$group, ref='NCFh0')
dds <- DESeq(dds)

#anno <- get_annotation(dds, 'hsapiens_gene_ensembl','ensembl_gene_id')
#anno <- na.omit(anno)
rld <- rlog(dds, blind=FALSE)

## TcGSA needs three things
## 1) kegg gene sets 
## 2) normalized data matrix
## 3) metadata

## Example code here: 
## https://github.com/borishejblum/TcGSA/issues/12

## Gene names in gene sets must match data matrix 
## 1 gene set

## C2 geneset downloaded from http://software.broadinstitute.org/gsea/msigdb/collections.jsp
TCGSA_C2_entrez <- GSA::GSA.read.gmt("c2.cp.v6.2.entrez.gmt")
TCGSA_kegg_entrez <- GSA::GSA.read.gmt("c2.cp.kegg.v6.2.entrez.gmt")


## 2 data matrix 
mat <- assay(rld)
genes <- row.names(mat)
genes.df <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

mat_tib <- 
  mat %>%
  as.tibble(rownames = "ENSEMBL") %>%
  left_join(genes.df, by = "ENSEMBL") %>%
  drop_na()

mat_tib_filt <- 
  mat_tib %>% 
  group_by(ENTREZID) %>% 
  filter(n()==1) %>% 
  ungroup() %>% 
  select(-ENSEMBL)

mat_mat <- as.matrix(mat_tib_filt[,1:48])
row.names(mat_mat) <- mat_tib_filt$ENTREZID

## 3 metadata

exp_design <- as.data.frame(colData(dds))
exp_design$Patient_ID <- as.factor(rep(c("P001","P001","P001","P001","P001","P001",
                       "P002","P002","P002","P002","P002","P002",
                       "P003","P003","P003","P003","P003","P003",
                       "P004","P004","P004","P004","P004","P004"),2))
exp_design$Sample_name <- row.names(design)
exp_design$timepoint <- as.numeric(rep(c('0','120','12','24','48','72'), 8))

##sanity checks
exp_design$Sample_name %in% colnames(mat_mat)
identical(ncol(mat_mat), nrow(exp_design))

## TcGSA analysis 

tcgsa_res <- TcGSA::TcGSA.LR(expr = mat_mat,
                             gmt = TCGSA_kegg_entrez,
                             design = exp_design,
                             subject_name = "Patient_ID",
                             time_name = "timepoint",
                             time_func = "cubic")
                             #group_name = "condition")



summary(tcgsa_res)

sig_genesets_adjPval <- as.tibble(TcGSA::signifLRT.TcGSA(tcgsa_res)$mixedLRTadjRes) %>% arrange(AdjPval) %>% filter(AdjPval < 0.05)

TcGSA::plot1GS(expr = as.data.frame(mat_mat),
     Subject_ID = exp_design$Patient_ID,
     TimePoint = exp_design$timepoint,
     gmt = TCGSA_kegg_entrez,
     geneset.name = "KEGG_PROTEASOME")


clust <- TcGSA::clustTrend(tcgs = tcgsa_res, 
                           expr=tcgsa_res$Estimations,
                           Subject_ID=exp_design$Patient_ID,
                           TimePoint=exp_design$timepoint,
                           baseline = 0,
                           group_of_interest="CF")
clust

plot(x=tcgsa_res, expr=as.data.frame(mat_mat),
     Subject_ID=exp_design$Patient_ID,
     TimePoint=exp_design$timepoint,
     group_of_interest="CF",
     clust_trends=clust,
     legend.breaks=seq(from=-2,to=2, by=0.01), time_unit="D",
     subtitle="test", cex.label.row=0.5, cex.label.col=1, cex.main=0.7,
     heatmap.width=0.2, dendrogram.size=0.3, margins=c(2,3),
     heatKey.size=0.8)
