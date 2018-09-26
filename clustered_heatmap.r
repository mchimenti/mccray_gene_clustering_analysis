## Create clustered heatmaps for timecourse and genes of interest
## Date: 09.18.2018
## Author: Michael Chimenti
## Organism: GRCh37 / human
## Aligners: salmon
## Design: 
## Reps: 4

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")
#biocLite("COMBINE-lab/wasabi")       #install wasabi the first time

#pseudoalignment analysis
library(wasabi)
library(sleuth)
library(biomaRt)
library(tidyverse)

setwd("~/iihg/RNA_seq/mccray/project_cf_rsv_sinnbrajesh_human_nov2017")

#########################
## Salmon > HDF5 > Sleuth
#########################

base_dir <- "~/iihg/RNA_seq/mccray/project_cf_rsv_sinnbrajesh_human_nov2017"

sal_dirs <- file.path(paste(base_dir, 
                            c("RSV_CF1_0hrs", "RSV_CF1_12hrs", "RSV_CF1_24hrs", "RSV_CF1_48hrs", "RSV_CF1_72hrs", "RSV_CF1_120hrs",
                              "RSV_CF2_0hrs", "RSV_CF2_12hrs", "RSV_CF2_24hrs", "RSV_CF2_48hrs", "RSV_CF2_72hrs", "RSV_CF2_120hrs",
                              "RSV_CF3_0hrs", "RSV_CF3_12hrs", "RSV_CF3_24hrs", "RSV_CF3_48hrs", "RSV_CF3_72hrs", "RSV_CF3_120hrs",
                              "RSV_CF4_0hrs", "RSV_CF4_12hrs", "RSV_CF4_24hrs", "RSV_CF4_48hrs", "RSV_CF4_72hrs", "RSV_CF4_120hrs",
                              "RSV_NCF1_0hrs", "RSV_NCF1_12hrs", "RSV_NCF1_24hrs", "RSV_NCF1_48hrs", "RSV_NCF1_72hrs", "RSV_NCF1_120hrs",
                              "MeV_NCF2_0hrs", "RSV_NCF2_12hrs", "RSV_NCF2_24hrs", "RSV_NCF2_48hrs", "RSV_NCF2_72hrs", "RSV_NCF2_120hrs",
                              "RSV_NCF3_0hrs", "RSV_NCF3_12hrs", "RSV_NCF3_24hrs", "RSV_NCF3_48hrs", "RSV_NCF3_72hrs", "RSV_NCF3_120hrs",
                              "RSV_NCF4_0hrs", "RSV_NCF4_12hrs", "RSV_NCF4_24hrs", "RSV_NCF4_48hrs", "RSV_NCF4_72hrs", "RSV_NCF4_120hrs"), "salmon", sep="/"))

prepare_fish_for_sleuth(sal_dirs)

## create a R matrix containing sample names and conditions from a text file
s2c <- read.table(file.path(base_dir, "design_sleuth.txt"), header = TRUE, stringsAsFactors=FALSE)

## add a column called "sal_dirs" containing paths to the data
s2c <- dplyr::mutate(s2c, path = sal_dirs)

## Get common gene names for transcripts

## this section queries Ensemble online database for gene names associated with transcript IDs
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- rename(t2g, 'target_id' = 'ensembl_transcript_id', 'ens_gene' = 'ensembl_gene_id', 'ext_gene' = 'external_gene_name')

## Create the Sleuth Object

## create the sleuth object with the information matrix and the transcript to gene mapping
so <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = "ens_gene", extra_bootstrap_summary = FALSE)

kaltab <- kallisto_table(so)
ifn_genes <- read.csv('~/iihg/RNA_seq/mccray/project_gene_clustering_cf_rsv_sept2018/typeI_interferon_inducible_genes.csv', header=FALSE)
names(ifn_genes) <- "ext_gene"
ifn_genes <- left_join(ifn_genes, t2g, on = 'ext_gene')
ifn_genes <- na.omit(ifn_genes) # drop not found gene names 
ifn_genes <- ifn_genes[ifn_genes$target_id %in% kaltab$target_id, ]  # drop transcripts not in the SO 

png("heatmap_IFN_transcripts_cluster.png", 1200, 1500, pointsize=20, res=100)
plot_transcript_heatmap(so, ifn_genes$target_id, units = 'tpm', cluster_rows = FALSE, 
                        color_low = "black", color_mid = "purple", color_high = "red",
                        show_rownames = FALSE, main = "~2800 IFN transcripts expr (log TPM), clustered by sample")
dev.off()

ifntab <- kaltab[kaltab$target_id %in% ifn_genes$target_id, ]

ifntab_mean <- ifntab %>% 
  as.tibble() %>%
  group_by(target_id) %>%
  mutate(mean_tpm = mean(tpm)) %>%
  arrange(desc(mean_tpm))
  
top60_ifn <- as.data.frame(unique(ifntab_mean[1:2800,]$target_id))
colnames(top60_ifn) <- "target_id"
top60_ifn <- left_join(top60_ifn, t2g, by = "target_id")

png("heatmap_IFN_transcripts_top60_cluster.png", 1200, 1500, pointsize=20, res=100)
plot_transcript_heatmap(so, top60_ifn$target_id, units = 'tpm', cluster_rows = TRUE, 
                        color_low = "black", color_mid = "purple", color_high = "red",
                        main = "Top 60 IFN transcripts by expr (log TPM), clustered by sample/exprs",
                        labels_row = top60_ifn$ext_gene)
dev.off()

####

# find DE transcripts between CF and nonCF 
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_tab_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_tab_tx <- dplyr::filter(sleuth_tab_tx, qval <= 0.01)

# subset only DE transcripts from IFN genes

ifn_DE <- ifn_genes[ifn_genes$target_id %in% sleuth_tab_tx$target_id,]

png("heatmap_IFN_transcripts_topDE_cluster.png", 1200, 1500, pointsize=20, res=100)
plot_transcript_heatmap(so, ifn_DE$target_id, units = 'tpm', cluster_rows = TRUE, 
                        color_low = "black", color_mid = "purple", color_high = "red",
                        main = "Top DE IFN transcripts by expr (log TPM), clustered by sample/exprs",
                        labels_row = ifn_DE$ext_gene)
dev.off()


