library(data.table)
library(parallel)
library(lavaan)
library(psych)
library(ggplot2)
library(introdataviz)
library(cowplot)

###### Phenotype prep
wgs_info.data <- fread("covariates/AMPAD_and_EADBBelgian_FullCovariatesAndSampleInfo_Correlation_19May2022.csv")

# Center gender and code female as 1, male as 0
wgs_info.data$gender <- as.numeric(as.factor(wgs_info.data$Sex)) - 1
wgs_info.data$gender <- 1-wgs_info.data$gender

# Specify continous age variable 
wgs_info.data$age <- wgs_info.data$AgeAtDeath_or_Sampling
wgs_info.data$age[wgs_info.data$age == "90_or_above"] <- 90
wgs_info.data$age <- as.numeric(wgs_info.data$age)

# Define ADvsCon as ordered for lavaan
wgs_info_pathological.data <- wgs_info.data[wgs_info.data$FK_final_CaseControl != "Other", ]
wgs_info_pathological.data$FK_final_CaseControl <- ordered(wgs_info_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Stratify by cohort
wgs_info_rosmap.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "ROSMAP_DLPFC", ]

# Load IDs for single nuclei RNAseq
id.data <- read.csv("data/rosmap_single_ID.csv")
wgs_info_rosmap_single.data <- merge(wgs_info_rosmap.data, id.data, by = "ID")

# Sinlge nuclei RNAseq data
snRNAseq.data <- read.table("data/ROSMAP_n465_NormalizedAvgExp_perCelltype_perIndividual.tsv")
snRNAseq.data <- as.data.frame(t(snRNAseq.data[1:4,]))

# Organize IDs and cell types in long format
snRNAseq.data$IID <- gsub('.*_','',row.names(snRNAseq.data))
snRNAseq.data$cell <- gsub('RNA.','',gsub('_.*','',row.names(snRNAseq.data)))

wgs_info_transcriptome_rosmap.data <- merge(wgs_info_rosmap_single.data, snRNAseq.data, by = "IID")

# Only expressed cell types
wgs_info_transcriptome_rosmap_expr_cell.data <- wgs_info_transcriptome_rosmap.data[
    wgs_info_transcriptome_rosmap.data$cell == "Ast" |
    wgs_info_transcriptome_rosmap.data$cell == "Exc" |
    wgs_info_transcriptome_rosmap.data$cell == "Inh" |
    wgs_info_transcriptome_rosmap.data$cell == "Mic" |
    wgs_info_transcriptome_rosmap.data$cell == "Oli" |
    wgs_info_transcriptome_rosmap.data$cell == "OPC",
]

### Plot distributions
## sex
EIF2S3 <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = EIF2S3, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("EIF2S3") +
  scale_fill_brewer(palette = "Accent")

KDM5C <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = KDM5C, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM5C") +
  scale_fill_brewer(palette = "Accent")

KDM6A <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = KDM6A, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM6A") +
  scale_fill_brewer(palette = "Accent")

PUDP <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = PUDP, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("PUDP") +
  scale_fill_brewer(palette = "Accent")

## AD
EIF2S3_ad <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = EIF2S3, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("EIF2S3") +
  scale_fill_brewer(palette = "Dark2")

KDM5C_ad <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM5C") +
  scale_fill_brewer(palette = "Dark2")

KDM6A_ad <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM6A") +
  scale_fill_brewer(palette = "Dark2")

PUDP_ad <- ggplot(wgs_info_transcriptome_rosmap_expr_cell.data, aes(x = cell, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("PUDP") +
  scale_fill_brewer(palette = "Dark2")

plots <- list(EIF2S3,EIF2S3_ad,KDM5C,KDM5C_ad,KDM6A,KDM6A_ad,PUDP,PUDP_ad)
all_genes_single <- plot_grid(plotlist = plots, ncol = 2)
ggsave(all_genes_single, file = "figures/all_genes_single.png", width = 5000, height = 5000, units = "px")


### SEM
mediation_rosmap.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER6 + ePEER7 + ePEER10 + ePEER18 + ePEER24 + ePEER40 + ePEER41 + ePEER42 + ePEER52 + ePEER60
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER6 + ePEER7 + ePEER10 + ePEER18 + ePEER24 + ePEER40 + ePEER41 + ePEER42 + ePEER52 + ePEER60
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

rosmap_transcripts <- c("EIF2S3","KDM5C","KDM6A","PUDP")

# Astrocytes
wgs_info_transcriptome_rosmap_ast.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Ast", ]

mediation_rosmap_ast.list <- mclapply(wgs_info_transcriptome_rosmap_ast.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_ast.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_ast.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_ast.data <- do.call(rbind, mediation_rosmap_ast.list)
mediation_rosmap_ast.data$cell <- "Astrocytes"

# Cell adhesion molecules
wgs_info_transcriptome_rosmap_cam.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "CAMs", ]

mediation_rosmap_cam.list <- mclapply(wgs_info_transcriptome_rosmap_cam.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_cam.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_cam.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_cam.data <- do.call(rbind, mediation_rosmap_cam.list)
mediation_rosmap_cam.data$cell <- "Cell Adhesion Molecules"

# Endothelial
wgs_info_transcriptome_rosmap_end.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "End", ]

mediation_rosmap_end.list <- mclapply(wgs_info_transcriptome_rosmap_end.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_end.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_end.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_end.data <- do.call(rbind, mediation_rosmap_end.list)
mediation_rosmap_end.data$cell <- "Endothelial"

# Excitatory neurons
wgs_info_transcriptome_rosmap_exc.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Exc", ]

mediation_rosmap_exc.list <- mclapply(wgs_info_transcriptome_rosmap_exc.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_exc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_exc.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_exc.data <- do.call(rbind, mediation_rosmap_exc.list)
mediation_rosmap_exc.data$cell <- "Excitatory Neurons"

# Fibroblast
wgs_info_transcriptome_rosmap_fib.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Fib", ]

mediation_rosmap_fib.list <- mclapply(wgs_info_transcriptome_rosmap_fib.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_fib.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_fib.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_fib.data <- do.call(rbind, mediation_rosmap_fib.list)
mediation_rosmap_fib.data$cell <- "Fibroblast"

# Inh
wgs_info_transcriptome_rosmap_inh.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Inh", ]

mediation_rosmap_inh.list <- mclapply(wgs_info_transcriptome_rosmap_inh.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_inh.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_inh.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_inh.data <- do.call(rbind, mediation_rosmap_inh.list)
mediation_rosmap_inh.data$cell <- "Inhibitory Neurons"

# Mic
wgs_info_transcriptome_rosmap_mic.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Mic", ]

mediation_rosmap_mic.list <- mclapply(wgs_info_transcriptome_rosmap_mic.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_mic.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_mic.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_mic.data <- do.call(rbind, mediation_rosmap_mic.list)
mediation_rosmap_mic.data$cell <- "Microglia"

# Oligodendrocyte
wgs_info_transcriptome_rosmap_oli.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Oli", ]

mediation_rosmap_oli.list <- mclapply(wgs_info_transcriptome_rosmap_oli.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_oli.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_oli.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_oli.data <- do.call(rbind, mediation_rosmap_oli.list)
mediation_rosmap_oli.data$cell <- "Oligodendrocytes"

# Oligodendrocyte precursor cells 
wgs_info_transcriptome_rosmap_opc.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "OPC", ]

mediation_rosmap_opc.list <- mclapply(wgs_info_transcriptome_rosmap_opc.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_opc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_opc.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_opc.data <- do.call(rbind, mediation_rosmap_opc.list)
mediation_rosmap_opc.data$cell <- "Oligodendrocyte Preursor Cells"

# Pericyte
wgs_info_transcriptome_rosmap_per.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "Per", ]

mediation_rosmap_per.list <- mclapply(wgs_info_transcriptome_rosmap_per.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_per.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_per.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_per.data <- do.call(rbind, mediation_rosmap_per.list)
mediation_rosmap_per.data$cell <- "Pericytes"

# Smooth Muscle Cells
wgs_info_transcriptome_rosmap_smc.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "SMC", ]

mediation_rosmap_smc.list <- mclapply(wgs_info_transcriptome_rosmap_smc.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_smc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_smc.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_smc.data <- do.call(rbind, mediation_rosmap_smc.list)
mediation_rosmap_smc.data$cell <- "Smooth Muscle Cells"

# T Cells
wgs_info_transcriptome_rosmap_T.data <- wgs_info_transcriptome_rosmap.data[wgs_info_transcriptome_rosmap.data$cell == "T", ]

mediation_rosmap_T.list <- mclapply(wgs_info_transcriptome_rosmap_T.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap_T.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap_T.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 4)
mediation_rosmap_T.data <- do.call(rbind, mediation_rosmap_T.list)
mediation_rosmap_T.data$cell <- "T Cells"

mediation_rosmap_snRNAseq.data <- rbind(mediation_rosmap_ast.data,mediation_rosmap_cam.data,mediation_rosmap_end.data,mediation_rosmap_exc.data,mediation_rosmap_fib.data,mediation_rosmap_inh.data,mediation_rosmap_mic.data,mediation_rosmap_oli.data,mediation_rosmap_opc.data,mediation_rosmap_per.data,mediation_rosmap_smc.data,mediation_rosmap_T.data)
write.csv(mediation_rosmap_snRNAseq.data, file = "results/mediation_rosmap_snRNAseq.csv")
