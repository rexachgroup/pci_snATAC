library(ArchR)
library(Matrix)
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/cisTopic"

projdir <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/projrmSubset_subPeak"
proj <- loadArchRProject(projdir)
dropSubCs <- c("undefined","ast.C6","mg.C8","mg.C10","mg.C15",
               "mg.C1", "mg.C5", "neu.C3", "neu.C4", "odc.C2", "neu.C2", "mg.C2")
proj$subClusters[proj$subClusters %like% "neuron"] <- gsub("neuron","neu",
                                                           proj$subClusters[proj$subClusters %like% "neuron"])
projclean <- proj[!proj$subClusters %in% dropSubCs,]


projclean$new_majorC <- projclean$subClusters
projclean$new_majorC[projclean$new_majorC %in% 
                       c("neu.C6","neu.C7","neu.C8","neu.C9")] <- "IN" 
projclean$new_majorC[projclean$new_majorC %like% "neu"] <- "EX" 
projclean$new_majorC <- gsub("\\..*","", projclean$new_majorC)
projclean <- projclean[!projclean$Sample %in% c("I1_7","P1_7_at1_7"),]
unique(projclean$Sample)
projclean #567510

# 
projmg <- projclean[projclean$subClusters %like% "mg",]

se <- getMatrixFromProject(projmg, useMatrix = "PeakMatrix")  # peaks x cells

mat <- assays(se)$PeakMatrix        # dgCMatrix
peak_gr <- rowRanges(se)            # GRanges of peaks
peak_names <- paste0(seqnames(peak_gr), ":", start(peak_gr), "-", end(peak_gr))
rownames(mat) <- peak_names

cell_md <- as.data.frame(colData(se))  # cell metadata

saveRDS(list(counts = mat, peaks = rownames(mat), cells = colnames(mat),
             meta = cell_md), file = file.path(outp, "Microglia_peaks_by_cells.rds"))
