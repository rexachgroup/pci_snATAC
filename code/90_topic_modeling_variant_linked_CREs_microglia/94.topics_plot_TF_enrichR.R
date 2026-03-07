# R 4.2.2
# XH Sep 9 2025
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(enrichR)
library(memes)
library(colorRamp2)
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/cisTopic"

peakmat <- readRDS(file.path(outp, "Microglia_peaks_by_cells.rds"))
meta <- peakmat$meta
dxs <- c("Control","AD","bvFTD","PSP_S")

# Pick Label
label <- "Dynamic"
#label <- "Dynamic_Stable"

topic_mat <- readRDS(file.path(outp, paste(label,"topic_matrices.rds",sep = ".")))

topic_by_cell <- topic_mat$topic_by_cell
topic_by_cre <- topic_mat$topic_by_cre

dim(topic_by_cre)

#"Assign" CREs to topics: top xx CREs per topic
N <- 70
top_cre_per_topic <- lapply(1:nrow(topic_by_cre), function(j) {
  ord <- order(topic_by_cre[j, ], decreasing = TRUE)
  colnames(topic_by_cre)[ord[1:N]]
})
names(top_cre_per_topic) <- rownames(topic_by_cre)



# 1 “Which CREs are in which topic?”
# z-score per CRE across topics not desired - instead scale per topic:
mat_scaled <- t(scale(t(topic_by_cre))) 

p <- Heatmap(mat_scaled,
             name = "NormTop (z)",
             cluster_rows = TRUE, cluster_columns = T,
             show_row_names = T,  show_column_names = F, 
             clustering_method_rows="ward.D2",
             clustering_method_column="ward.D2",
             column_title = "Topic contribution per CRE",
             col = colorRamp2(c(-2,0,2), c("#4575b4","#f7f7f7","#d73027")))

pdf(file.path(outp, paste("heat.topics_x_CREs",label,"pdf",sep = ".")), height = 4)
p
dev.off()

#check the NormTop per topic
for (i in names(top_cre_per_topic)){
  cres <- top_cre_per_topic[[i]]
  m <- mat_scaled[i, cres]
  print(min(m)) # min > 0.19
}

# Save to file
write.csv(stack(top_cre_per_topic),
          file = file.path(outp, paste(label,"Assigned_CREs_perTopic.csv",sep = ".")),
          row.names = FALSE)

# 2 “Which Topics contribute to disease/subc/dx_x_subc group”
plot_group_heatmap <- function(avg_dx){
  mat <- as.matrix(avg_dx)
  # scale rows to highlight relative topic usage per group
  mat_row_scaled <- scale(t(mat))
  p <- pheatmap(mat_row_scaled,
                clustering_method = "ward.D2", scale = "none",
                col = colorRamp2(c(-2,0,2), c("#2c67b0","#f7f7f7","#b4222d")),
                main = "Topics scaled contribution")
  return(p)
}

# dx
avg_dx <- NULL
for (dx in c("Control","AD","bvFTD","PSP_S")){
  cells <- meta[meta$Clinical.Dx %in% dx,] %>% rownames()
  cells_keep <- intersect(cells, colnames(topic_by_cell))
  mat_ss <- topic_by_cell[,cells_keep]
  avg_dx <- cbind(rowMeans(mat_ss), avg_dx)
}
colnames(avg_dx) <- c("Control","AD","bvFTD","PSP_S")

# subc
avg_subc <- NULL
subCs <- unique(meta$subClusters) %>% sort()
for (subc in subCs){
  cells <- meta[meta$subClusters %in% subc,] %>% rownames()
  cells_keep <- intersect(cells, colnames(topic_by_cell))
  mat_ss <- topic_by_cell[,cells_keep]
  avg_subc <- cbind(rowMeans(mat_ss), avg_subc)
}
colnames(avg_subc) <- subCs

# subc-dxs
avg_subcDx <- NULL
subCs <- unique(meta$subClusters) %>% sort()
for (subc in subCs){
  for (dx in dxs){
    cells <- meta[meta$subClusters %in% subc &
                    meta$Clinical.Dx %in% dx,] %>% rownames()
    cells_keep <- intersect(cells, colnames(topic_by_cell))
    mat_ss <- topic_by_cell[,cells_keep]
    avg_subcDx <- cbind(rowMeans(mat_ss), avg_subcDx)
  }
}
colnames(avg_subcDx) <- apply(expand.grid(subCs, dxs), 1, paste0, collapse = "-")

#
pdf(file.path(outp, paste(label,"heat.topics_x_Dx.pdf",sep = ".")), height = 3, width = 5)
plot_group_heatmap(avg_dx)
dev.off()

pdf(file.path(outp, paste(label,"heat.topics_x_subCs.pdf",sep = ".")), height = 4, width = 5)
plot_group_heatmap(avg_subc)
dev.off()

pdf(file.path(outp, paste(label,"heat.topics_x_subCdx.pdf",sep = ".")), height = 8, width = 5)
plot_group_heatmap(avg_subcDx)
dev.off()

pdf(file.path(outp, paste(label,"heat.topics_x_subCdx_each.pdf",sep = ".")), height = 3, width = 5)
for (subc in subCs){
  avg_ss <- avg_subcDx[,colnames(avg_subcDx) %like% subc]
  plot_group_heatmap(avg_ss) %>% print()
}
dev.off()


# 3. TFs and enrichment for topics
# load CRE annotation "pciATAC_peakset"
source("/geschwindlabshares/RexachGroup/Xia_Data/My_R_libs/pci_snATAC_peakSet_w_Enh.R")
cre_gene <- pciATAC_peakset[,c("peakType","nearestGene")] %>% unique()
cre_gene$peak <- paste(cre_gene@seqnames, ":", cre_gene@ranges@start, 
                        "-", cre_gene@ranges@start+500, sep = "")

topic_genes <- list()
gr_topic_cre <- list()
for (topic in names(top_cre_per_topic)){
  cres <- top_cre_per_topic[[topic]]
  genes <- cre_gene[cre_gene$peak %in% cres,]$nearestGene %>% unique()
  topic_genes[[topic]] <- genes
  
  gr <- GRanges(seqnames = gsub("^([^:]+):.*", "\\1", cres),
                ranges   = IRanges(start = as.numeric(gsub(".*\\:(.*)-.*", "\\1", cres)),
                                   end   = as.numeric(gsub(".*-(.*)", "\\1", cres) ))
  )
  gr$peak <- cres
  
  m <- findMatches(gr, cre_gene)

  gr$peakType <- rep(NA, length(gr))
  gr$peakGene <-  rep(NA, length(gr))
  gr[queryHits(m)]$peakType <- as.character(pciATAC_peakset[subjectHits(m)]$peakType)
  gr[queryHits(m)]$peakGene <- as.character(pciATAC_peakset[subjectHits(m)]$nearestGene)

  gr_topic_cre[[topic]] <- gr
  
  # check
  #setdiff(genes,gr$peakGene)
  #setdiff(gr$peakGene,genes)

}



# 3.1 enrichR
library(enrichR)
dbs <- listEnrichrDbs()
dbs$libraryName
usedbs <- c("GO_Biological_Process_2021","KEGG_2021_Human")
#         "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")

all.result <- NULL #all save in one df
for(group in names(topic_genes)){
  genes <- topic_genes[[group]]
  print(paste(group, ":", length(genes)))
  
  enriched <- enrichr(genes, usedbs)
  
  for (i in 1:length(usedbs)){
    db <- usedbs[i]
    result <- enriched[[db]]
    # add two cols: "group" ""db"
    annote_group <- rep(group, dim(result)[1])
    annote_db <- rep(db, dim(result)[1])
    result <- cbind(group=annote_group, db=annote_db, result)
    all.result <- rbind(all.result, result)
  }
}

write.csv(all.result, file.path(outp, paste(label,"enrichR_topics.csv",sep = ".")))

# plot selected enriched terms
# picked_terms (from excel) saved in 'picked_enrichedR_topics.txt'
enriched_term <- read.csv(file.path(outp, "Dynamic.enrichR_topics.csv"))
dt_term <- read.csv(file.path(outp, "picked_enrichedR_topics.txt"),sep = "\t")
pick_term <- unique(dt_term$Term)

dropterms <- c("positive regulation of cysteine-type endopeptidase activity involved in apoptotic process",
               "neutrophil activation involved in immune response",
               "positive regulation of intracellular signal transduction",
               "antigen receptor−mediated signaling pathway",
               "positive regulation of phagocytosis",
               "positive regulation of cytokine production involved in immune response"
               )
dtp_enrich <- enriched_term[enriched_term$Term %in% pick_term, ]
dtp_enrich <- dtp_enrich[dtp_enrich$Adjusted.P.value < 0.1, ]
dtp_enrich$nGenes <- gsub("\\/.*","",dtp_enrich$Overlap) %>% as.numeric()
gsub(" \\(.*","",dtp_enrich$Term)
dtp_enrich$Term <- gsub(" \\(.*","",dtp_enrich$Term)

intersect(dropterms, dtp_enrich$Term)
dtp_enrich <- dtp_enrich[!dtp_enrich$Term %in% dropterms, ]

dtp_enrich$Term <- factor(dtp_enrich$Term, levels=unique(dtp_enrich$Term))

p_term <- ggplot(dtp_enrich, aes(x = group, y = Term)) +
  geom_point(aes(colour = -log10(Adjusted.P.value), size = nGenes)) +
  labs(x="",y="", title = "Enriched terms for genes linking to CREs per topic")+
  #scale_colour_gradientn(colours=ArchR::ArchRPalettes$solarExtra) +
  scale_colour_gradient(low = "pink", high = "darkred") +
  scale_size(range=c(2,5)) +
  theme_set(theme_bw())+
  theme(plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,face = "plain"),
        axis.title.y=element_text(size=10,vjust = 2,hjust = 0.5,face = "plain"),
        axis.text.x=element_text(size=10, face="plain", vjust = 1,hjust = 1, angle = 45),
        axis.text.y=element_text(size=10,face="plain"),
        legend.position = "right",
        #legend.position = c(1,1),
        #legend.title=element_blank(),
        legend.text=element_text(size=10),
        panel.grid =element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"))

pdf(file.path(outp, paste(label,"enrichedTerms_topics.pdf",sep = ".")),
    width = 9, height = 6)
p_term
dev.off()


# 3.2 meems - TF enrichment
library(memes)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
OUTBASE <- "/home/xhan/data_cellrangerOut/project0318/outdir_projrmSubset_subPeak/CRE_projrmSubset_subPeak/motif_by_memes"
PATHDB<- "/geschwindlabshares/RexachGroup/SharedData/motif_databases"
#motifdb <- file.path(PATHDB, "JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme")
#db_use <- "JASPAR" # only JASPAR  not for CIS-BP HUMAN
motifdb <- file.path(PATHDB, "CIS-BP_2.00/Homo_sapiens.meme")
db_use <- "CIS-BP_2.00"
memeBin <- "/home/xhan/apps/meme-5.4.1/bin"

query_seq <- get_sequence(gr_topic_cre, genome)
dir <- file.path(outp, "TF_Topic_DynCRE-frVar-QTL")
if (!dir.exists(dir)){
  dir.create(dir)
}

#run enrichment
for (group in names(gr_topic_cre)){
  print(group)
  dir_tf <- file.path(dir, group)
  if (!dir.exists(dir)){
    dir.create(dir)
  }
  enrichMotif <- runAme(input = query_seq[[group]],
                        control = "shuffle", #InvariantSeqs[[group]]
                        outdir = dir_tf,
                        method = "fisher", # as both input and control are 501bp
                        database = motifdb,
                        meme_path = memeBin)
}

#--
fetch_memesOut <- function(group, dir_n){
  df <- read.csv(file.path(dir_n, "ame.tsv"), sep = "\t")
  df <- df[!df$motif_DB %in% "",]
  
  # drop ZNF* TFs and found in other species
  df <- df[!df$motif_alt_ID %like% "ZNF",]
  df <- df[!df$motif_alt_ID %like% "\\(",]
  df <- df[complete.cases(df),]
  #groupn <- str_match(group, 'P_(.*)?_unique$')[,2]
  ameOut <- data.frame(group = rep(group, dim(df)[1]),
                       motif = df[["motif_alt_ID"]],
                       Padj_mlog10 = -log10(df[["adj_p.value"]]))
  ameOut$scaled <- scale(ameOut$Padj_mlog10)
  return(ameOut)
}

plot_p <- function(dtp, titlen){
  # use only the scaled_mPadj > 0 & make other tfs NA as 0
  #dtp <- dtp[dtp$scaled > 0, ]
  order_motifs <- unique(dtp$motif)
  
  dtp$motif <- factor(dtp$motif, levels = order_motifs)
  add <-  theme_set(theme_bw())+
    theme( plot.title = element_text(size=12, face = "bold", hjust=0.5),
           axis.title.x=element_text(size=12,face = "plain"),
           title =element_text(size=12, face='bold'),
           axis.title.y=element_text(size=12,vjust = 2,hjust = 0.5,face = "plain"),
           axis.text.x=element_text(size=12, face="plain", angle = 90, vjust = 0.5),
           #axis.text.x=element_blank(),
           axis.text.y=element_text(size=12,face="plain"),
           legend.position = "right",
           #legend.title=element_blank(),
           legend.text=element_text(size=9),
           panel.grid = element_line(linewidth = 0.1),
           #panel.grid.minor.y = element_line(linewidth=1, colour = "grey"),
           #panel.border = element_rect(linetype = "dashed", fill = NA),
           #panel.border = element_rect(color = "black", fill = NA, size = 0.1),
           panel.border = element_blank(),
           axis.ticks = element_line(colour = "black", linewidth= 0.1),
           #axis.line = element_line(linewidth= 0.1, colour = "black", linetype="solid")
           axis.line = element_blank())
  my_palette <- colorRampPalette(ArchR::ArchRPalettes$whitePurple)(100)
  p <- ggplot(dtp, aes(x = group, y = motif, fill = scaled)) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colours=my_palette) +
    #coord_flip() +
    #labs(x="",y="",fill = "-log10(P-adj)", title = titlen) + add
    labs(x="",y="",fill = "Norm -log10(P-adj)", title = titlen) + add
  p2 <- ggplot(dtp, aes(x = group, y = motif, fill = Padj_mlog10)) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colours=my_palette) +
    #coord_flip() +
    #labs(x="",y="",fill = "-log10(P-adj)", title = titlen) + add
    labs(x="",y="",fill = "     -log10(P-adj)", title = titlen) + add
  return(list(p,p2))
}

dt_ameOut <- NULL
for (group in rownames(topic_by_cre)){
  dir_tf <- file.path(dir, group)
  dt_ameOut <- rbind(dt_ameOut, fetch_memesOut(group, dir_tf))
}

pdf(file.path(outp, paste(label,"enrichedTF_topics.pdf",sep = ".")),
    width = 5, height = 15)
plot_p(dt_ameOut, "TFs enriched in Topics")
dev.off()

# plot mg-related TFs 
microglia_TFs <- c(
  "RFX5",   # MHC-II antigen presentation
  "JUNB", "FOSL2", "CREB1", "MYC", "MAX", "EGR1", # AP-1 / stress / proliferation
  "HSF1", "TFE3", "YY1", # stress, lysosome, chromatin looping
  "IRF1", "IRF7", "SPI1", "SPIB", "SPIC", # innate immune & interferon
  "CEBPA", "CEBPB", "CEBPD", "CEBPG", "CEBPE", # inflammatory CEBPs
  "KLF2", "KLF4", "KLF6", "KLF7", "KLF9", "KLF10", "KLF12", "KLF13", "KLF15", # immune/metabolic regulation
  "MEF2A", "MEF2B", "MEF2C", "MEF2D", # homeostatic suppression of activation
  "BACH2", "RARA", # immune regulation, retinoic acid signaling
  "SP1", "SP4", # general immune/chromatin regulation
  "CLOCK", "ARNT2", # circadian & hypoxia responses
  "ETV6", "ELF1", "ELF2", "ELF5", "FLI1", # ETS family, immune-related
  "BCL11A", # immune development & chromatin regulation
  "ERG1", # topic 2
  "USF1","VSX2", "ZFP69B", "RHOXF2", "HKR1" # add topic 9
)

drop_TFs <- c("ELF1","KLF15","ELF5","ELF1","KLF13","KLF7","YY1","HSF1","CEBPD","KLF12")
microglia_TFs<-microglia_TFs[!microglia_TFs %in% drop_TFs]

pdf(file.path(outp, paste(label,"enrichedTF_topics_picked.pdf",sep = ".")),
    width = 5, height = 6.5)
plot_p(dt_ameOut[dt_ameOut$motif %in% microglia_TFs,], "(Microglia-related) TFs enriched in Topics")
dev.off()

# 3.3 TFBS for each CRE
# annotate TFBS for peaks given TFs - by MEME's FIMO
for (group in names(gr_topic_cre)){
  fa_file <- paste("DynCRE",group,"fa",sep = ".")
  writeXStringSet(query_seq[[1]], filepath = file.path(outp, "fimo_output",fa_file))
}


# Run the FIMO to find TFBS (for all TFs)
setwd(file.path(outp, "fimo_output"))
for (group in names(gr_topic_cre)){
  print(group)
  fa_file <- paste("DynCRE",group,"fa",sep = ".")
  outdir <- group
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  fimo_command <- paste(
    paste("/geschwindlabshares/RexachGroup/Xia_Data/tool/bin/bin/fimo --oc",outdir, motifdb, fa_file)
  )
  
  system(fimo_command)
}
