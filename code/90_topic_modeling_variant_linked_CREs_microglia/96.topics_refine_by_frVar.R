# Summarized topics with SNPs-3 types

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
CREs_w_SNPs <- read.csv(file.path(outp, "Mg_CREs_w_frVar_eQTL_wide.csv"), row.names = 1)
dynCRE_w_SNPs <- CREs_w_SNPs[CREs_w_SNPs$peakset %in% "Dynamic",]

# 2. Make data frame
df_topics <- NULL
for (topic in topics){
  v_cres <- topic_cre[topic_cre$ind %in% topic, ]$values
  dt_cres_w_snp <- CREs_w_SNPs[CREs_w_SNPs$loci %in% v_cres,]
  dt_cres_w_snp$topic <- topic
  df_topics <- rbind(df_topics, dt_cres_w_snp)
}

#-----write csv, add TFBS for the cre
dt <- df_topics %>%
  mutate(
    topic      = as.character(topic),
    gene       = as.character(peakGene),
    peak_id    = as.character(loci)      # CRE coordinates
  )

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

write.csv(combined, file.path(outp, "topics_w_snps_TF.csv"))
