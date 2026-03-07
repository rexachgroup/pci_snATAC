#conda activate cistopic
# ***!!!! mutate .bashrc R 4.2.2 !!! all versions of R

# Identify modules for CREs contain frVar+eQTL

# XH Sep 8 2025
library(data.table)
library(dplyr)
library(GenomicRanges)
library(cisTopic)

outp <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/cisTopic"
dxs <- c("Control","AD","bvFTD","PSP_S")

# 1. Load peakmat and CREs
# Peak mat of microglia (from 01.*.R) : 900k peaks x 50k cells
peakmat <- readRDS(file.path(outp, "Microglia_peaks_by_cells.rds"))
mat <- peakmat$counts
meta <- peakmat$meta

# Microglia (specific) CREs overlapped with frVars and QTL (from 02.*.R)
CREs_w_SNPs <- read.csv(file.path(outp, "Mg_CREs_w_frVar_eQTL.csv"), row.names = 1)
table(CREs_w_SNPs$peakset)

dynCRE <- CREs_w_SNPs[CREs_w_SNPs$peakset %in% "Dynamic",]$loci
stableCRE <- CREs_w_SNPs[CREs_w_SNPs$peakset %in% "Stable",]$loci




# Pick Label
label <- "Dynamic"
#label <- "Dynamic_Stable"




if(label %in% "Dynamic"){
  keep <- rownames(mat) %in% dynCRE
}else{
  keep <- rownames(mat) %in% c(dynCRE, stableCRE)
}

mat_keep <- mat[keep, , drop = FALSE]

# Binarize input (cisTopic accepts binary)
mat_keep@x[mat_keep@x > 0] <- 1

# 2. Run cisTopic (can try several times to pick a good K based on figure!!)
cisTopicObject <- createcisTopicObject(count.matrix = mat_keep,    # or mat_bin
                                       project.name   = label,
                                       min.cells      = 1,
                                       min.regions    = 1,
                                       keepCountsMatrix = TRUE)

cisTopicObject <- runWarpLDAModels(cisTopicObject,
                                   topic=c(2:30, 35, 40),
                                   nCores     = 8,
                                   seed       = 123,
                                   iterations = 1000)


pdf(file.path(outp, paste(label,"select_models.pdf", sep = ".")))
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()


# select K and save
cisTopicObject <- selectModel(cisTopicObject, type='derivative') 

topic_by_cell <- modelMatSelection(cisTopicObject, 
                                   target = "cell",  method = "Probability")
topic_by_cre <- modelMatSelection(cisTopicObject, 
                                  target = "region",  method = "NormTop", all.regions = TRUE)
topic_mat <- list(
  "topic_by_cell"=topic_by_cell,
  "topic_by_cre"=topic_by_cre
)

save(cisTopicObject, file=file.path(outp, paste0(label,".cisTopicObject_full.RData")))
saveRDS(topic_mat, file.path(outp, paste(label,"topic_matrices.rds",sep = ".")))
                                  