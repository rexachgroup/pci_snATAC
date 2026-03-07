# CREs with MPRA frVar & sneQTL
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ArchR)
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/cisTopic"

# 1.1 ATAC peak set
# load CRE annotation
source("/geschwindlabshares/RexachGroup/Xia_Data/My_R_libs/pci_snATAC_peakSet_w_Enh.R")

# load stable and dynamic peaks of microglia, from (ldsc input)
path_mypeak <- "/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/mypeaks"
bed_files <- c("mg_dxMarkerP_ctUniqueP.bed", "mg_dxInvariantP_ctUniqueP.bed")
list_atacpeaks <- list()
for (file in bed_files){
  fname <- ifelse(file %like% "InvariantP", "Stable","Dynamic")
  dt <- read.table(file.path(path_mypeak, file))
  gr <- GRanges(seqnames = dt$V1, ranges = IRanges(dt$V2, end = dt$V3))
  
  # add PE-gene info
  m <- findMatches(gr, pciATAC_peakset)
  gr$peakType <- rep(NA, length(gr))
  gr$peakGene <-  rep(NA, length(gr))
  gr[queryHits(m)]$peakType <- as.character(pciATAC_peakset[subjectHits(m)]$peakType)
  gr[queryHits(m)]$peakGene <- as.character(pciATAC_peakset[subjectHits(m)]$nearestGene)
  
  list_atacpeaks[[fname]] <- gr
}
names(list_atacpeaks)


# 2.1 Load two eQTL datasets
eQTL_Malhotra_w_pos <- readRDS("/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/out/eQTL_Malhotra_wpos.rds")
eQTL_Malhotra_mg <- eQTL_Malhotra_w_pos[eQTL_Malhotra_w_pos$celltype %in% "Microglia",]
gr_qtl_Malhotra <- GRanges(seqnames = paste("chr",eQTL_Malhotra_mg$gene_chr, sep = ""), 
                           ranges = IRanges(eQTL_Malhotra_mg$snp_pos_hg38, 
                                            end = eQTL_Malhotra_mg$snp_pos_hg38), 
                           snp=eQTL_Malhotra_mg$snp, 
                           QTLgene = eQTL_Malhotra_mg$eGene,
                           beta=eQTL_Malhotra_mg$beta,
                           pval=eQTL_Malhotra_mg$pval)


eQTL_DeJager_mg <- NULL
ct <- "Mic"
inp_1 <- "/geschwindlabshares/RexachGroup/SharedData/QTL/sneQT_Fujita_and_DeJager_2024/celltype"
filen <- paste("celltype-eqtl-sumstats",ct,"tsv.gz" ,sep = ".")
qtl <- read.csv(file.path(inp_1, filen), sep = "\t")
sig_qtl <- qtl[qtl$significant_by_2step_FDR %in% "Yes",]
eQTL_DeJager_mg <- rbind(eQTL_DeJager_mg, sig_qtl)

gr_qtl_DeJager <- GRanges(seqnames = eQTL_DeJager_mg$chr38, 
                           ranges = IRanges(eQTL_DeJager_mg$pos38, 
                                            end = eQTL_DeJager_mg$pos38), 
                          snp=eQTL_DeJager_mg$snps, 
                          QTLgene = eQTL_DeJager_mg$gene_symbol,
                          beta=eQTL_DeJager_mg$beta,
                          pval=eQTL_DeJager_mg$pvalue)

# 2.2 Load MPRA frVar
mpra <- readRDS("/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/out/obj.mpra_w_LDblock.rds")
colnames(mpra)[4] <- "mpra_variant_id"
mpra <- mpra[!is.na(mpra$start_hg38), ]

used_snp <-  mpra[mpra$pval < 0.05,]
gr_frVar <- GRanges(
  seqnames=paste("chr", used_snp$chr_hg38, sep = ""),
  IRanges(start=used_snp$start_hg38,
          end=used_snp$start_hg38+1)
)
values(gr_frVar) <- used_snp

# 3. Get CREs with SNPs
count_Peak_wQTL <- function(gr_snp, gr_peak){
    o <- findOverlaps(gr_peak, gr_snp)
    
    tmp_peak <- as.data.frame(gr_peak)
    gr_peak$loci <- paste(tmp_peak$seqnames, ":", tmp_peak$start, "-", tmp_peak$end, sep = "")
    
    df.overlap <- data.frame(gr_peak[queryHits(o)], gr_snp[subjectHits(o)])
    df.overlap <- df.overlap[,c("loci","peakType","peakGene","snp","QTLgene","pval","beta")]
    df.overlap <- df.overlap[!is.na(df.overlap$peakGene), ]
    #  count overlap only when eGene=peakGene
    keep.overlap <- df.overlap[df.overlap$peakGene == df.overlap$QTLgene,]
    keep.overlap$QTLgene <- NULL
  return(keep.overlap)
}

count_Peak_wMPRA <- function(gr_snp, gr_peak){
  o <- findOverlaps(gr_peak, gr_snp)
  
  tmp_peak <- as.data.frame(gr_peak)
  gr_peak$loci <- paste(tmp_peak$seqnames, ":",tmp_peak$start, "-",tmp_peak$end, sep = "")
  
  df.overlap <- data.frame(gr_peak[queryHits(o)], gr_snp[subjectHits(o)])
  df.overlap <- df.overlap[,c("loci","peakType","peakGene","mpra_variant_id",
                              "pval","logFC")]
  colnames(df.overlap)[4] <- "snp"
  colnames(df.overlap)[6] <- "beta"
  return(df.overlap)
}

dt_peaks_w_SNPs <- NULL
for (peakset in names(list_atacpeaks)){
  gr_peak <- list_atacpeaks[[peakset]]
  dt1 <- count_Peak_wQTL(gr_qtl_Malhotra, gr_peak)
  dt1$snp_type <- "eQTL_Malhotra"
  dt2 <- count_Peak_wQTL(gr_qtl_DeJager, gr_peak)
  dt2$snp_type <- "eQTL_DeJager"
  dt3 <- count_Peak_wMPRA(gr_frVar, gr_peak)
  dt3$snp_type <- "MPRA_frVar"
  dt <- rbind(dt1, dt2, dt3)
  dt$peakset <- peakset
  dt_peaks_w_SNPs <- rbind(dt_peaks_w_SNPs, dt)
}

CREs_w_SNPs <- dt_peaks_w_SNPs[dt_peaks_w_SNPs$peakType %in% c("PE","Promoter"), ]
table(CREs_w_SNPs$peakset)

write.csv(CREs_w_SNPs, file.path(outp, "Mg_CREs_w_frVar_eQTL.csv"))
