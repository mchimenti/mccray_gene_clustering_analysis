library("TcGSA")
library("GEOquery")

temp <- tempfile()
utils::download.file("http://doi.org/10.1371/journal.pcbi.1004310.s007", destfile = temp, mode = "wb")
load(unz(temp, "ReproducibleRFiles/GMTs_PLOScb.RData", open = "r"))
unlink(temp)
rm(temp)

GEOquery::getGEOSuppFiles('GSE46734')

design <- read.delim(gzfile("GSE46734/GSE46734_DALIA1longitudinalTranscriptome_DESIGN_anonym.txt.gz"))
design_preATI <- design[-which(design$TimePoint<0 | design$TimePoint==16 | design$TimePoint>22), ]
head(design_preATI,5)

expr_preATI <- read.delim(gzfile("GSE46734/GSE46734_DALIA1longitudinalTranscriptome_PALO01_PreATI_NEQC_NonParamCombat.txt.gz"))
rownames(expr_preATI) <- expr_preATI$PROBE_ID
expr_preATI <- expr_preATI[,as.character(design_preATI$Sample_name)]

expr_preATI[1:4,1:4]

identical(ncol(expr_preATI), nrow(design_preATI))

tcgsa_result <- TcGSA::TcGSA.LR(expr = expr_preATI, 
                                gmt = gmt_modulesV2, 
                                design = design_preATI, 
                                subject_name = "Patient_ID", 
                                time_name = "TimePoint")
summary(tcgsa_result)

TcGSA::plot1GS(expr = expr_preATI, 
               #plot1GS(expr = tcgsa_result$Estimations,
               gmt = gmt_modulesV2, 
               Subject_ID = design_preATI$Patient_ID, 
               TimePoint = design_preATI$TimePoint,
               clustering = FALSE, 
               time_unit = "W", 
               geneset.name = "M3.2", 
               title="",
               margins=0.4,
               lab.cex=0.37,
               axis.cex=0.37,
               line.size=0.45,
               gg.add=list(ggplot2::theme(legend.position="none"),
                           ggplot2::ylim(-1.26,1.26)
               ))

clust <- TcGSA::clustTrend(tcgs = tcgsa_result, 
                           expr=tcgsa_result$Estimations,
                           Subject_ID=design_gen$Patient_ID,
                           TimePoint=design_gen$TimePoint,
                           baseline = 0)
clust

plot(x=tcgsa_result, expr=tcgsa_result$Estimations,
     Subject_ID=design_preATI$Patient_ID,
     TimePoint=design_preATI$TimePoint,
     clust_trends=clust,
     legend.breaks=seq(from=-2,to=2, by=0.01), time_unit="D",
     subtitle="Pneumo vs Saline", cex.label.row=0.5, cex.label.col=1, cex.main=0.7,
     heatmap.width=0.2, dendrogram.size=0.3, margins=c(2,3),
     heatKey.size=0.8)
