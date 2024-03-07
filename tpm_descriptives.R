library(data.table)
library(ggplot2)
library(introdataviz)
library(cowplot)

# Hits to describe
hits <- c("ENSG00000005889.15",
          "ENSG00000126012.11",
          "ENSG00000130021.13",
          "ENSG00000130741.10",
          "ENSG00000147050.14",
          "ENSG00000169895.5",
          "ENSG00000198034.10",
          "ENSG00000224975.1",
          "ENSG00000229807.9",
          "ENSG00000270223.1")

# ROSMAP
tpm_rosmap.data <- fread("../amp_ad_genotypes/gct/ROSMAP_n560_rnaseqc_stranded_tpm_06June2020.gct.gz", data.table = F)
tpm_rosmap.data <- tpm_rosmap.data[tpm_rosmap.data$Name %in% hits, ]
rownames(tpm_rosmap.data) <- tpm_rosmap.data$Description
tpm_rosmap_t.data <- as.data.frame(t(tpm_rosmap.data[,3:562]))
tpm_rosmap_t.data$Region <- "DLPFC"
tpm_rosmap_id.data <- fread("../amp_ad_genotypes/gct/ROSMAP_n560_sample_participation_lookup.tsv", data.table = F)
tpm_rosmap_t.data <- merge(tpm_rosmap_t.data, tpm_rosmap_id.data, by.x = "row.names", by.y = "sample_id")

# MAYO
tpm_mayo.data <- fread("../amp_ad_genotypes/gct/MayoRNASeq_TCX_n259_rnaseqc_tpm_12january2020.gct.gz", data.table = F)
tpm_mayo.data <- tpm_mayo.data[tpm_mayo.data$Name %in% hits, ]
rownames(tpm_mayo.data) <- tpm_mayo.data$Description
tpm_mayo_t.data <- as.data.frame(t(tpm_mayo.data[,3:261]))
tpm_mayo_t.data$Region <- "TCX"
tpm_mayo_id.data <- fread("../amp_ad_genotypes/gct/MayoRNASeq_TCX_n259_sample_participation_lookup.tsv", data.table = F)
tpm_mayo_t.data <- merge(tpm_mayo_t.data, tpm_mayo_id.data, by.x = "row.names", by.y = "sample_id")

# MSBB BA10
tpm_ba10.data <- fread("../amp_ad_genotypes/gct/MSBB_BA10_n207_rnaseqc_tpm_12january2020.gct.gz", data.table = F)
tpm_ba10.data <- tpm_ba10.data[tpm_ba10.data$Name %in% hits, ]
rownames(tpm_ba10.data) <- tpm_ba10.data$Description
tpm_ba10_t.data <- as.data.frame(t(tpm_ba10.data[,3:209]))
tpm_ba10_t.data$Region <- "BA10"
tpm_ba10_id.data <- fread("../amp_ad_genotypes/gct/MSBB_BA10_n207_sample_participation_lookup.tsv", data.table = F)
tpm_ba10_t.data <- merge(tpm_ba10_t.data, tpm_ba10_id.data, by.x = "row.names", by.y = "sample_id")

# MSBB BA22
tpm_ba22.data <- fread("../amp_ad_genotypes/gct/MSBB_BA22_n186_rnaseqc_tpm_12january2020.gct.gz", data.table = F)
tpm_ba22.data <- tpm_ba22.data[tpm_ba22.data$Name %in% hits, ]
rownames(tpm_ba22.data) <- tpm_ba22.data$Description
tpm_ba22_t.data <- as.data.frame(t(tpm_ba22.data[,3:188]))
tpm_ba22_t.data$Region <- "BA22"
tpm_ba22_id.data <- fread("../amp_ad_genotypes/gct/MSBB_BA22_n186_sample_participation_lookup.tsv", data.table = F)
tpm_ba22_t.data <- merge(tpm_ba22_t.data, tpm_ba22_id.data, by.x = "row.names", by.y = "sample_id")

# MSBB BA36
tpm_ba36.data <- fread("../amp_ad_genotypes/gct/MSBB_BA36_n162_rnaseqc_tpm_12january2020.gct.gz", data.table = F)
tpm_ba36.data <- tpm_ba36.data[tpm_ba36.data$Name %in% hits, ]
rownames(tpm_ba36.data) <- tpm_ba36.data$Description
tpm_ba36_t.data <- as.data.frame(t(tpm_ba36.data[,3:164]))
tpm_ba36_t.data$Region <- "BA36"
tpm_ba36_id.data <- fread("../amp_ad_genotypes/gct/MSBB_BA36_n162_sample_participation_lookup.tsv", data.table = F)
tpm_ba36_t.data <- merge(tpm_ba36_t.data, tpm_ba36_id.data, by.x = "row.names", by.y = "sample_id")

# MSBB BA10
tpm_ba44.data <- fread("../amp_ad_genotypes/gct/MSBB_BA44_n178_rnaseqc_tpm_12january2020.gct.gz", data.table = F)
tpm_ba44.data <- tpm_ba44.data[tpm_ba44.data$Name %in% hits, ]
rownames(tpm_ba44.data) <- tpm_ba44.data$Description
tpm_ba44_t.data <- as.data.frame(t(tpm_ba44.data[,3:180]))
tpm_ba44_t.data$Region <- "BA44"
tpm_ba44_id.data <- fread("../amp_ad_genotypes/gct/MSBB_BA44_n178_sample_participation_lookup.tsv", data.table = F)
tpm_ba44_t.data <- merge(tpm_ba44_t.data, tpm_ba44_id.data, by.x = "row.names", by.y = "sample_id")

# Merge all regions
tpm.data <- rbind(tpm_rosmap_t.data,tpm_mayo_t.data,tpm_ba10_t.data,tpm_ba22_t.data,tpm_ba36_t.data,tpm_ba44_t.data)

# Add clinical info
wgs_info.data <- fread("covariates/AMPAD_and_EADBBelgian_FullCovariatesAndSampleInfo_Correlation_19May2022.csv", data.table = F)
wgs_info_pathological.data <- wgs_info.data[wgs_info.data$FK_final_CaseControl != "Other", ]
wgs_info_pathological.data$FK_final_CaseControl <- ordered(wgs_info_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Remove duplicate information on MSBB samples, due to repeated measures
wgs_info_pathological.data <- wgs_info_pathological.data[!duplicated(wgs_info_pathological.data$ID), ]

tpm.data <- merge(tpm.data, wgs_info_pathological.data[c("ID","Sex","FK_final_CaseControl")], by.x = "participation_id", by.y = "ID")

tpm.data$region_sex <- as.factor(paste(tpm.data$Region,tpm.data$Sex))
tpm.data$Region <- factor(tpm.data$Region, levels = c("DLPFC","TCX","BA10","BA22","BA36","BA44"))
levels(tpm.data$Region) <- c("DLPFC (ROS/MAP)","TCX (Mayo)","BA10 (MSBB)","BA22 (MSBB)","BA36 (MSBB)","BA44 (MSBB)")

### Violin plot
## Sex-stratified
# EIF2S3
EIF2S3 <- ggplot(tpm.data, aes(x = Region, y = EIF2S3, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("EIF2S3")

# KDM5C
KDM5C <- ggplot(tpm.data, aes(x = Region, y = KDM5C, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("KDM5C")

# RPS4X
RPS4X <- ggplot(tpm.data, aes(x = Region, y = RPS4X, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("RPS4X")

# KDM6A
KDM6A <- ggplot(tpm.data, aes(x = Region, y = KDM6A, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  ggtitle("KDM6A") + scale_y_sqrt()

# PUDP
PUDP <- ggplot(tpm.data, aes(x = Region, y = PUDP, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("PUDP")

# ZFX
ZFX <- ggplot(tpm.data, aes(x = Region, y = ZFX, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  ggtitle("ZFX")

# RP13-36G14.4
RP13 <- ggplot(tpm.data, aes(x = Region, y = `RP13-36G14.4`, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  ggtitle("RP13-36G14.4")

# SYAP1
SYAP1 <- ggplot(tpm.data, aes(x = Region, y = SYAP1, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("SYAP1")

# XIST
XIST <- ggplot(tpm.data, aes(x = Region, y = XIST, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (log10 transformed)") + scale_y_continuous(trans='log10') +
  ggtitle("XIST")

# INE1
INE1 <- ggplot(tpm.data, aes(x = Region, y = INE1, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Accent") +
  ggtitle("INE1")
  scale_fill_brewer(palette = "Dark2")

## AD status
# EIF2S3
EIF2S3_ad <- ggplot(tpm.data, aes(x = Region, y = EIF2S3, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("EIF2S3")

# KDM5C
KDM5C_ad <- ggplot(tpm.data, aes(x = Region, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM5C")

# RPS4X
RPS4X_ad <- ggplot(tpm.data, aes(x = Region, y = RPS4X, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RPS4X")

# KDM6A
KDM6A_ad <- ggplot(tpm.data, aes(x = Region, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM6A") + scale_y_sqrt()

# PUDP
PUDP_ad <- ggplot(tpm.data, aes(x = Region, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  ggtitle("PUDP")

# ZFX
ZFX_ad <- ggplot(tpm.data, aes(x = Region, y = ZFX, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  ggtitle("ZFX")

# RP13-36G14.4
RP13_ad <- ggplot(tpm.data, aes(x = Region, y = `RP13-36G14.4`, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RP13-36G14.4")

# SYAP1
SYAP1_ad <- ggplot(tpm.data, aes(x = Region, y = SYAP1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("SYAP1")

# XIST
XIST_ad <- ggplot(tpm.data, aes(x = Region, y = XIST, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (log10 transformed)") + scale_y_continuous(trans='log10') +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("XIST") +
  guides(fill=guide_legend(title="AD status")) +
  scale_fill_brewer(palette = "Dark2")

# INE1
INE1_ad <- ggplot(tpm.data, aes(x = Region, y = INE1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  ggtitle("INE1") +
  guides(fill=guide_legend(title="AD status")) +
  scale_fill_brewer(palette = "Dark2")
  
plots <- list(EIF2S3,EIF2S3_ad,KDM5C,KDM5C_ad,RPS4X,RPS4X_ad,KDM6A,KDM6A_ad,PUDP,PUDP_ad,ZFX,ZFX_ad,RP13,RP13_ad,SYAP1,SYAP1_ad,XIST,XIST_ad,INE1,INE1_ad)
all_genes <- plot_grid(plotlist = plots, ncol = 2)
ggsave(all_genes, file = "figures/all_genes.png", width = 4000, height = 10000, units = "px")

