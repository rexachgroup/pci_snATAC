# Summarized topics with SNPs info (to one table or visualization)

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(colorRamp2)
library(ggplotify)
library(patchwork)
library(gridExtra)
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/cisTopic"

add <-   theme_set(theme_bw())+
  theme( plot.title = element_text(size=12, face = "bold"),
         plot.subtitle = element_text(size = 10, face = "italic", color = "black", hjust = 0.5),
         axis.title.x=element_text(size=12,face = "plain"),
         title =element_text(size=12, face='bold'),
         axis.title.y=element_text(size=12,vjust = 2,hjust = 0.5,face = "plain"),
         axis.text.x=element_text(size=12, face="plain", angle = 90, hjust = 1, color = "black"),
         #axis.text.x=element_blank(),
         axis.text.y=element_text(size=12,face="plain", color = "black"),
         legend.position = "right",
         legend.title=element_text(size=12),
         legend.text=element_text(size=12),
         panel.grid = element_blank(),
         #panel.grid.minor.y = element_line(linewidth=1, colour = "grey"),
         #panel.border = element_rect(linetype = "dashed", fill = NA),
         panel.border = element_blank(),
         axis.ticks = element_line(colour = "black", linewidth= 0.1),
         axis.line = element_line(linewidth= 0.1, colour = "black", linetype="solid"))


# 1. Load data
# cres per topic
topic_cre <- read.csv(file.path(outp, "Dynamic.Assigned_CREs_perTopic.csv"))
topics <- unique(topic_cre$ind)

colors.topics <- ArchR::ArchRPalettes$circus[1:length(topics)]
names(colors.topics) <- topics

# tfbs of cres per topic (MEF2C JUNB)
tfbs <- read.csv(file.path(outp, "fimo_output","DynCRE_Topics_TFBS.csv"),row.names = 1)
tfbs[tfbs$peak %in% peak,]

# frVar eQTL in CREs
CREs_w_SNPs <- read.csv(file.path(outp, "Mg_CREs_w_frVar_eQTL.csv"), row.names = 1)
dynCRE_w_SNPs <- CREs_w_SNPs[CREs_w_SNPs$peakset %in% "Dynamic",]

# FTD GWAS pval
path_gwas <- "/geschwindlabshares/RexachGroup/Xia_Data/resources/GWAS/FTD_META_bvFTDpval_updated_rsID_w_hg38_loci.txt"
FTD_gwas <- read.csv(path_gwas, sep = ",", row.names = 1)
FTD_gwas[!is.na(FTD_gwas$SNP), ] %>% dim() # 6,020,472
FTD_gwas <- FTD_gwas[!is.na(FTD_gwas$SNP), ]
dim(FTD_gwas)

#rownames(FTD_gwas) <- FTD_gwas$SNP

# 2. Make data frame
df_topics <- NULL
for (topic in topics){
  v_cres <- topic_cre[topic_cre$ind %in% topic, ]$values
  dt_cres_w_snp <- CREs_w_SNPs[CREs_w_SNPs$loci %in% v_cres,]
  # add gwas pval and beta
  with_gwas <- FTD_gwas[FTD_gwas$SNP %in% dt_cres_w_snp$snp, ]
  rownames(with_gwas) <- with_gwas$SNP
  with_gwas <- with_gwas[dt_cres_w_snp$snp, ]
  dt_cres_w_snp$gwas_pval <- with_gwas$pValue
  dt_cres_w_snp$gwas_pval_bvFTD <- with_gwas$pval_bvFTDgwas
  dt_cres_w_snp$gwas_beta <- with_gwas$beta1
  dt_cres_w_snp$topic <- topic
  df_topics <- rbind(df_topics, dt_cres_w_snp)
}


dt_show <- df_topics[!is.na(df_topics$gwas_pval), ]

#---3. SNP distribution & Plot 
table(dt_show$topic)

# (1)  Distribution of SNPs pvalue
# Functional Pval/Beta and even GWAS-Pval of SNPs across Topics
# x - topics, y - gwas_pval, boxplot dot plot

p1 <- ggplot(dt_show, aes(x = topic, y = -log10(gwas_pval), colour = topic)) +
  geom_jitter(size=1) +
  #geom_boxplot(outlier.shape = 21, alpha = 0.6) +
  scale_colour_manual(values=colors.topics) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dotted") +
  labs(x = "", y = "-log10(eQTL/MPRA p-value)",
       title = "Distribution of P values for frVar and eQTL in each topic") +
  add
p2 <- ggplot(dt_show, aes(x = topic, y = -log10(pval), colour = topic)) +
  geom_jitter(size=1) +
  #geom_boxplot(outlier.shape = 21, alpha = 0.6) +
  scale_colour_manual(values=colors.topics) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dotted") +
  labs(x = "", y = "-log10(FTD GWAS p-value)",
       title = "Distribution of GWAS-pvalues for frVar and eQTL in each topic") +
  add
p3 <- ggplot(dt_show, aes(x = topic, y = beta, colour = topic)) +
  geom_jitter(size=1) +
  #geom_boxplot(outlier.shape = 21, alpha = 0.6) +
  scale_colour_manual(values=colors.topics) +
  labs(x = "", y = "beta (eQTL/MPRA)",
       title = "Distribution of GWAS p values for frVar and eQTL in each topic") +
  add

pdf(file.path(outp, paste("pDot.topics_SNPs.DynCRE.FTD_gwas.pdf",sep = ".")), height = 4)
p1
p2
p3
dev.off()

# Plot (2) Per topic, functional variants P value vs GWAS pval 
ps <- list()
for (topic in topics){
  dtp <- dt_show[dt_show$topic %in% topic, ]
  dtp <- dtp[,c("snp","peakGene","gwas_pval","pval","topic")]
  print(paste(topic, "#sig snps:", length(unique(dtp$snp)), 
              "#sig genes", length(unique(dtp$peakGene)) ))
  
  p_scatter <- ggplot(dtp, aes(x = -log10(pval), y = -log10(gwas_pval), colour = topic)) +
    geom_jitter(size=1) +
    #geom_boxplot(outlier.shape = 21, alpha = 0.6) +
    scale_colour_manual(values=colors.topics) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dotted") +
    geom_vline(xintercept = -log10(5e-8), color = "red", linetype = "dotted") +
    # add labels only for SNPs passing gwas-p cutoff
    ggrepel::geom_text_repel(
      data = subset(dtp, pval < 5e-8 | gwas_pval < 5e-8),  # only significant SNPs
      aes(label = peakGene),   
      size = 2,
      max.overlaps = 20,   # allow more labels
      box.padding = 0.5,   # extra space around text
      point.padding = 0.2, # space between text and point
      segment.color = "grey50", # line color
      segment.size = 0.4,       # line thickness
      segment.alpha = 0.8       # line transparency
    ) +
    labs(x = "-log10(eQTL/MPRA p-value)", y = "-log10(FTD GWAS p-value)",
         title = paste("SNP Distribution in", topic) )+
    add
  
  # table of significant genes (sorted by pval)
  df_tab <- dtp %>%
    filter(!is.na(pval), pval < 5e-8 | gwas_pval < 5e-8) %>%
    arrange(pval) %>%
    transmute(Gene = peakGene,
              `-log10(p)` = sprintf("%.2f", -log10(pval)))
  
  if (nrow(df_tab) == 0) {
    tab_plot <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "No SNPs pass cutoff",
               size = 1) +
      labs(title = "")
  } else {
    tab_grob  <- tableGrob(df_tab, rows = NULL,
                           theme = ttheme_default(
                             core = list(
                               fg_params = list(cex = 0.6, lineheight = 0.8)  # text smaller & tighter
                             ),
                             colhead = list(
                               fg_params = list(cex = 0.65, fontface = "bold") # header smaller
                             ),
                             padding = unit(c(1, 1), "mm")  # reduce vertical/horizontal cell padding
                           )
    )
    tab_plot  <- as.ggplot(tab_grob) + labs(title = "")
  }
  
  # combine: scatter | table
  ps[[topic]] <- (p_scatter | tab_plot) + plot_layout(widths = c(1,1))
}

pdf(file.path(outp, paste("pDot.PxP_per_topic.DynCRE.FTD_gwas.pdf",sep = ".")), height = 3, width = 6)
ps
dev.off()

# bvFTD pval
ps_bvftd <- list()
for (topic in topics){
  dtp <- dt_show[dt_show$topic %in% topic, ]
  dtp <- dtp[,c("snp","peakGene","gwas_pval_bvFTD","pval","topic")]
  print(paste(topic, "#sig snps:", length(unique(dtp$snp)), 
              "#sig genes", length(unique(dtp$peakGene)) ))
  
  p_scatter <- ggplot(dtp, aes(x = -log10(pval), y = -log10(gwas_pval_bvFTD), colour = topic)) +
    geom_jitter(size=1) +
    #geom_boxplot(outlier.shape = 21, alpha = 0.6) +
    scale_colour_manual(values=colors.topics) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dotted") +
    geom_vline(xintercept = -log10(5e-8), color = "red", linetype = "dotted") +
    # add labels only for SNPs passing gwas-p cutoff
    ggrepel::geom_text_repel(
      data = subset(dtp, pval < 5e-8 | gwas_pval_bvFTD < 5e-8),  # only significant SNPs
      aes(label = peakGene),   
      size = 2,
      max.overlaps = 20,   # allow more labels
      box.padding = 0.5,   # extra space around text
      point.padding = 0.2, # space between text and point
      segment.color = "grey50", # line color
      segment.size = 0.4,       # line thickness
      segment.alpha = 0.8       # line transparency
    ) +
    labs(x = "-log10(eQTL/MPRA p-value)", y = "-log10(FTD GWAS p-value)",
         title = paste("SNP Distribution in", topic) )+
    add
  
  # table of significant genes (sorted by pval)
  df_tab <- dtp %>%
    filter(!is.na(pval), pval < 5e-8 | gwas_pval_bvFTD < 5e-8) %>%
    arrange(pval) %>%
    transmute(Gene = peakGene,
              `-log10(p)` = sprintf("%.2f", -log10(pval)))
  
  if (nrow(df_tab) == 0) {
    tab_plot <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "No SNPs pass cutoff",
               size = 1) +
      labs(title = "")
  } else {
    tab_grob  <- tableGrob(df_tab, rows = NULL,
                           theme = ttheme_default(
                             core = list(
                               fg_params = list(cex = 0.6, lineheight = 0.8)  # text smaller & tighter
                             ),
                             colhead = list(
                               fg_params = list(cex = 0.65, fontface = "bold") # header smaller
                             ),
                             padding = unit(c(1, 1), "mm")  # reduce vertical/horizontal cell padding
                           )
    )
    tab_plot  <- as.ggplot(tab_grob) + labs(title = "")
  }
  
  # combine: scatter | table
  ps_bvftd[[topic]] <- (p_scatter | tab_plot) + plot_layout(widths = c(1,1))
}

pdf(file.path(outp, paste("pDot.PxP_per_topic.DynCRE.bvFTD_gwas.pdf",sep = ".")), height = 3, width = 6)
ps_bvftd
dev.off()

#-----write csv, add TFBS for the cre
# plot gene - cre (x n) - tf (x n) - snps (x n)
dt_show <- df_topics[!is.na(df_topics$gwas_pval), ]
dt_show$mlog10_pval <- -log10(dt_show$pval)
dt_show$mlog10_gwasP <- -log10(dt_show$gwas_pval)

dt <- dt_show %>%
  mutate(
    topic      = as.character(topic),
    gene       = as.character(peakGene),
    peak_id    = as.character(loci)      # CRE coordinates
  ) %>%
  select(topic, gene, peak_id, peakType, snp, pval, beta, snp_type, peakset, 
         gwas_pval, gwas_beta, mlog10_pval, mlog10_gwasP)

tf <- tfbs %>%
  mutate(
    topic      = as.character(group),
    gene       = as.character(peakGene),
    peak_id    = as.character(peak)      # should match dt$loci
  ) %>%
  # Keep only one row per TFBS hit (de-dup in case of overlaps within the same peak)
  distinct(topic, gene, peak_id, TF, motif_match, .keep_all = TRUE) %>%
  select(topic, gene, peak_id, TF, motif_match)

intersect(dt$peak_id, tf$peak_id)

combined <- merge(dt, tf, by = c("topic", "peak_id"), all.x = TRUE)

#check
peak <- "chr4:155674005-155674505"
tf[tf$peak_id %in% peak,]
dt[dt$peak_id %in% peak,]
combined[combined$peak_id %in% peak, ]

write.csv(combined, file.path(outp, "table_topics_snps_TF.FTD_gwas.csv"))


# 4. Fisher test for snps enriched in topics

df <- df_topics[!is.na(df_topics$gwas_pval), ]

# Mark “genome-wide significant” 
alpha <- 5e-8
df$sig <- FALSE
#df[df$pval < alpha | df$gwas_pval < alpha,]$sig <- TRUE # for both gwas-p and funtional p
df[df$pval < alpha,]$sig <- TRUE # for only functional p
df[,c("pval","gwas_pval","sig")]

# Totals over ALL topics (background)
N_total <- n_distinct(df$snp)
N_sig   <- n_distinct(df$snp[df$sig])

# Per-topic counts and Fisher’s test
dt_fisher <- NULL
for (topic in topics){
  n_topic <- df[df$topic %in% topic,]$snp %>% unique() %>% length()
  n_sig_topic <- df[df$topic %in% topic & df$sig %in% TRUE,]$snp %>% unique() %>% length()
  n_notsig_nottopic <- df[!df$topic %in% topic & 
                            df$sig %in% FALSE,]$snp %>% unique() %>% length()
  
  # 2×2 table: in topic vs not, significant vs not
  mat <- matrix(c(
    n_sig_topic,                      n_topic - n_sig_topic,
    N_sig - n_sig_topic,              n_notsig_nottopic
  ), byrow = T, nrow = 2)
  
  rownames(mat) <- c("in_topic", "not_in_topic")
  colnames(mat) <- c("sig_snp","not_sig_snp")
  print(mat)
  
  fisher_out <- fisher.test(mat)
  p_fisher <- fisher_out$p.value
  or  <- fisher_out$estimate
  # the odds ratio (OR) from Fisher’s exact test is basically a relative measure of how 
  #enriched a topic is for significant SNPs compared to the rest of the genome/background.
  
  dto <- data.frame(topic, p_fisher, OR=or)
  dt_fisher <- rbind(dt_fisher, dto)
}

dt_fisher$fdr <- p.adjust(dt_fisher$p_fisher, method = "BH")
dt_fisher


pv <- ggplot(dt_fisher, aes(x = log2(OR), y = -log10(fdr), label = topic)) +
  geom_point(aes(color = (fdr < 0.1))) +
  scale_color_manual(values = c("TRUE"="red","FALSE"="grey40")) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggrepel::geom_text_repel(size = 4) +
  labs(title = "Topics enriched in SNPs by Fisher test",
       x = "log2 Odds ratio (enrichment)",
       y = "-log10(FDR)") +
  add
pdf(file.path(outp, paste("pDot.fisher_topic_w_snps.DynCRE.FTD_gwas.pdf",sep = ".")), height = 3, width = 4)
pv
dev.off()
