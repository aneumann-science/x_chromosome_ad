library(data.table)
library(ggplot2)
library(introdataviz)
library(cowplot)

# Hits to describe
hits <- c("ENSG00000126012.11",
          "ENSG00000130021.13",
          "ENSG00000147050.14",
          "ENSG00000169895.5",
          "ENSG00000198034.10",
          "ENSG00000225470.6",
          "ENSG00000234969.1")

# ROSMAP
tpm_rosmap.data <- fread("expression_matrices_and_PEER/TPMs/ROSMAP_n636_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_rosmap.data <- tpm_rosmap.data[tpm_rosmap.data$Name %in% hits, ]
rownames(tpm_rosmap.data) <- tpm_rosmap.data$Description
tpm_rosmap_t.data <- as.data.frame(t(tpm_rosmap.data[,3:638]))
tpm_rosmap_t.data$Region <- "DLPFC"
row.names(tpm_rosmap_t.data) <- gsub("Aligned","",row.names(tpm_rosmap_t.data))
row.names(tpm_rosmap_t.data) <- gsub("redo","",row.names(tpm_rosmap_t.data))

# MAYO
tpm_mayo.data <- fread("expression_matrices_and_PEER/TPMs/MayoRNASeq_n264_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_mayo.data <- tpm_mayo.data[tpm_mayo.data$Name %in% hits, ]
rownames(tpm_mayo.data) <- tpm_mayo.data$Description
tpm_mayo_t.data <- as.data.frame(t(tpm_mayo.data[,3:266]))
tpm_mayo_t.data$Region <- "TCX"
row.names(tpm_mayo_t.data) <- gsub("Aligned","",row.names(tpm_mayo_t.data))

# MSBB BA10
tpm_ba10.data <- fread("expression_matrices_and_PEER/TPMs/MSBB_BA10_n227_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_ba10.data <- tpm_ba10.data[tpm_ba10.data$Name %in% hits, ]
rownames(tpm_ba10.data) <- tpm_ba10.data$Description
tpm_ba10_t.data <- as.data.frame(t(tpm_ba10.data[,3:229]))
tpm_ba10_t.data$Region <- "BA10"
row.names(tpm_ba10_t.data) <- gsub("Aligned","",row.names(tpm_ba10_t.data))

# MSBB BA22
tpm_ba22.data <- fread("expression_matrices_and_PEER/TPMs/MSBB_BA22_n202_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_ba22.data <- tpm_ba22.data[tpm_ba22.data$Name %in% hits, ]
rownames(tpm_ba22.data) <- tpm_ba22.data$Description
tpm_ba22_t.data <- as.data.frame(t(tpm_ba22.data[,3:204]))
tpm_ba22_t.data$Region <- "BA22"
row.names(tpm_ba22_t.data) <- gsub("Aligned","",row.names(tpm_ba22_t.data))

# MSBB BA36
tpm_ba36.data <- fread("expression_matrices_and_PEER/TPMs/MSBB_BA36_n180_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_ba36.data <- tpm_ba36.data[tpm_ba36.data$Name %in% hits, ]
rownames(tpm_ba36.data) <- tpm_ba36.data$Description
tpm_ba36_t.data <- as.data.frame(t(tpm_ba36.data[,3:182]))
tpm_ba36_t.data$Region <- "BA36"
row.names(tpm_ba36_t.data) <- gsub("Aligned","",row.names(tpm_ba36_t.data))

# MSBB BA44
tpm_ba44.data <- fread("expression_matrices_and_PEER/TPMs/MSBB_BA44_n193_rnaseqc_stranded_tpm_14April2025.gct.gz", data.table = F)
tpm_ba44.data <- tpm_ba44.data[tpm_ba44.data$Name %in% hits, ]
rownames(tpm_ba44.data) <- tpm_ba44.data$Description
tpm_ba44_t.data <- as.data.frame(t(tpm_ba44.data[,3:195]))
tpm_ba44_t.data$Region <- "BA44"
row.names(tpm_ba44_t.data) <- gsub("Aligned","",row.names(tpm_ba44_t.data))

# Merge all regions
tpm.data <- rbind(tpm_rosmap_t.data,tpm_mayo_t.data,tpm_ba10_t.data,tpm_ba22_t.data,tpm_ba36_t.data,tpm_ba44_t.data)

# Add clinical info
pheno.data <- fread("phenotype/amp_ad_pheno.csv")

# Define ADvsCon as ordered for lavaan
pheno_pathological.data <- pheno.data[pheno.data$FK_final_CaseControl != "Other", ]
pheno_pathological.data$FK_final_CaseControl <- ordered(pheno_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Remove duplicate information on MSBB samples, due to repeated measures
pheno_pathological.data <- pheno_pathological.data[!duplicated(pheno_pathological.data$specimenID), ]

tpm.data <- merge(tpm.data, pheno_pathological.data[,c("specimenID","Sex","FK_final_CaseControl")], by.x = "row.names", by.y = "specimenID")

tpm.data$Sex <- factor(tpm.data$Sex, levels = c("Male","Female"))
tpm.data$Region <- factor(tpm.data$Region, levels = c("DLPFC","TCX","BA10","BA22","BA36","BA44"))
levels(tpm.data$Region) <- c("DLPFC (ROS/MAP)","TCX (Mayo)","BA10 (MSBB)","BA22 (MSBB)","BA36 (MSBB)","BA44 (MSBB)")

### Violin plot
## Sex-stratified
# KDM5C
KDM5C <- ggplot(tpm.data, aes(x = Region, y = KDM5C, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") +
  ggtitle("KDM5C")

# PUDP
PUDP <- ggplot(tpm.data, aes(x = Region, y = PUDP, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") +
  ggtitle("PUDP")

# KDM6A
KDM6A <- ggplot(tpm.data, aes(x = Region, y = KDM6A, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  ggtitle("KDM6A") + scale_y_sqrt()

# SYAP1
SYAP1 <- ggplot(tpm.data, aes(x = Region, y = SYAP1, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  ggtitle("SYAP1") + scale_y_sqrt()

# RPS4X
RPS4X <- ggplot(tpm.data, aes(x = Region, y = RPS4X, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") +
  ggtitle("RPS4X")

# JPX
JPX <- ggplot(tpm.data, aes(x = Region, y = JPX, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  ggtitle("JPX")

# RP13-33H18.1
RP13 <- ggplot(tpm.data, aes(x = Region, y = `RP13-33H18.1`, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  ggtitle("RP13-33H18.1")

## AD status
# KDM5C
KDM5C_ad <- ggplot(tpm.data, aes(x = Region, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM5C")

# PUDP
PUDP_ad <- ggplot(tpm.data, aes(x = Region, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  ggtitle("PUDP")

# KDM6A
KDM6A_ad <- ggplot(tpm.data, aes(x = Region, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)")  + scale_fill_brewer(palette = "Dark2") + scale_y_sqrt() +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM6A")

# SYAP1
SYAP1_ad <- ggplot(tpm.data, aes(x = Region, y = SYAP1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("SYAP1") + scale_y_sqrt()

# RPS4X
RPS4X_ad <- ggplot(tpm.data, aes(x = Region, y = RPS4X, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RPS4X")

# JPX
JPX_ad <- ggplot(tpm.data, aes(x = Region, y = JPX, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("JPX") + scale_y_sqrt()

# RP13-33H18.1
RP13_ad <- ggplot(tpm.data, aes(x = Region, y = `RP13-33H18.1`, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RP13-33H18.1") + scale_y_sqrt()
  
plots <- list(KDM5C,KDM5C_ad,PUDP,PUDP_ad,KDM6A,KDM6A_ad,SYAP1,SYAP1_ad,RPS4X,RPS4X_ad,JPX,JPX_ad,RP13,RP13_ad)
all_genes <- plot_grid(plotlist = plots, ncol = 2)
ggsave(all_genes, file = "figures/all_genes.png", width = 4000, height = 10000, units = "px")

###### Plots by sex
## Male
# KDM5C
KDM5C_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM5C Male")

# PUDP
PUDP_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("TPM") +
  ggtitle("PUDP Male")

# KDM6A
KDM6A_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_y_sqrt() +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM6A Male")

# SYAP1
SYAP1_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = SYAP1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("SYAP1 Male") + scale_y_sqrt()

# RPS4X
RPS4X_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = RPS4X, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RPS4X Male")

# JPX
JPX_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = JPX, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("JPX Male") + scale_y_sqrt()

# RP13-33H18.1
RP13_male <- ggplot(tpm.data[tpm.data$Sex == "Male",], aes(x = Region, y = `RP13-33H18.1`, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RP13-33H18.1 Male") + scale_y_sqrt()

## Female
# KDM5C
KDM5C_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM5C Female")

# PUDP
PUDP_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  ggtitle("PUDP Female")

# KDM6A
KDM6A_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)")  + scale_fill_brewer(palette = "Dark2") + scale_y_sqrt() +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("KDM6A Female")

# SYAP1
SYAP1_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = SYAP1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("SYAP1 Female") + scale_y_sqrt()

# RPS4X
RPS4X_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = RPS4X, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RPS4X Female")

# JPX
JPX_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = JPX, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("JPX Female") + scale_y_sqrt()

# RP13-33H18.1
RP13_female <- ggplot(tpm.data[tpm.data$Sex == "Female",], aes(x = Region, y = `RP13-33H18.1`, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("TPM (square root transformed)") + scale_fill_brewer(palette = "Dark2") +
  guides(fill=guide_legend(title="AD status")) +
  ggtitle("RP13-33H18.1 Female") + scale_y_sqrt()

plots_sex <- list(KDM5C_male,KDM5C_female,PUDP_male,PUDP_female,KDM6A_male,KDM6A_female,SYAP1_male,SYAP1_female,RPS4X_male,RPS4X_female,JPX_male,JPX_female,RP13_male,RP13_female)
all_genes_sex <- plot_grid(plotlist = plots_sex, ncol = 2)
ggsave(all_genes_sex, file = "figures/all_genes_sex.png", width = 4000, height = 10000, units = "px")

