library(data.table)
library(parallel)
library(lavaan)
library(psych)
library(ggplot2)
library(introdataviz)
library(cowplot)

###### Phenotype prep
pheno.data <- fread("phenotype/amp_ad_pheno.csv")

# Center gender and code female as 1, male as 0
pheno.data$gender <- as.numeric(as.factor(pheno.data$Sex)) - 1
pheno.data$gender <- 1-pheno.data$gender

# Specify continuous age variable 
pheno.data$age <- pheno.data$AgeAtDeath
pheno.data$age[pheno.data$age == "90_or_above"] <- 90
pheno.data$age <- as.numeric(pheno.data$age)

# Define ADvsCon as ordered for lavaan
pheno_pathological.data <- pheno.data[pheno.data$FK_final_CaseControl != "Other", ]
pheno_pathological.data$FK_final_CaseControl <- ordered(pheno_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Load imputed genetic components
load("phenotype/gPC_imputed.Rdata")
pheno_pathological.data <- merge(pheno_pathological.data, gPC_imputed.data, by = "specimenID")

# Stratify by cohort
pheno_rosmap.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "ROSMAP_DLPFC", ]

## Load peer factors
rosmap_peer.data <- fread("expression_matrices_and_PEER/ePEER_60_ROSMAP_DLPFC.autosomes.allSamples_n636.14April2025.PEER_covariates.txt", data.table = F)

# Prepare for merging
peer_prep <- function(peer.data) {
  rownames(peer.data) <- peer.data$ID
  peer.data <- peer.data[,2:ncol(peer.data)]
  peer.data <- as.data.frame(t(peer.data))
  peer.data$specimenID <- rownames(peer.data)
  return(peer.data)
}

rosmap_peer.data <- peer_prep(rosmap_peer.data)

pheno_rosmap.data <- merge(pheno_rosmap.data, rosmap_peer.data, by = "specimenID")

# Sinlge nuclei RNAseq data
snRNAseq.data <- fread("expression_matrices_and_PEER/ROSMAP_n465_NormalizedAvgExp_perCelltype_perIndividual_5May2025.tsv")
snRNAseq.data <- data.frame(snRNAseq.data, row.names=1)
snRNAseq.data <- as.data.frame(t(snRNAseq.data[1:5,]))

# Organize IDs and cell types in long format
snRNAseq.data$IID <- gsub('.*_','',row.names(snRNAseq.data))
snRNAseq.data$cell <- gsub('RNA.','',gsub('_.*','',row.names(snRNAseq.data)))

# Merge with phenotype data
pheno_snc_rosmap.data <- merge(pheno_rosmap.data, snRNAseq.data, by.x = "individualID", by.y = "IID")
# Number of participants
dim(unique(pheno_snc_rosmap.data[pheno_snc_rosmap.data$FK_final_CaseControl == "AD", "individualID"]))
dim(unique(pheno_snc_rosmap.data[pheno_snc_rosmap.data$FK_final_CaseControl == "Control", "individualID"]))

# Only expressed cell types
pheno_snc_rosmap.data_expr_cell.data <- pheno_snc_rosmap.data[
    pheno_snc_rosmap.data$cell == "Ast" |
    pheno_snc_rosmap.data$cell == "Exc" |
    pheno_snc_rosmap.data$cell == "Inh" |
    pheno_snc_rosmap.data$cell == "Mic" |
    pheno_snc_rosmap.data$cell == "Oli" |
    pheno_snc_rosmap.data$cell == "OPC",
]

### Plot distributions
## sex
KDM5C <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = KDM5C, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM5C") +
  scale_fill_brewer(palette = "Accent")

PUDP <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = PUDP, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("PUDP") +
  scale_fill_brewer(palette = "Accent")

KDM6A <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = KDM6A, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM6A") +
  scale_fill_brewer(palette = "Accent")

SYAP1 <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = SYAP1, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("SYAP1") +
  scale_fill_brewer(palette = "Accent")

JPX <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = JPX, fill = Sex)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  ylab("Log Normalized Expression") +
  ggtitle("JPX") +
  scale_fill_brewer(palette = "Accent")

## AD
KDM5C_ad <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = KDM5C, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM5C") +
  scale_fill_brewer(palette = "Dark2")

PUDP_ad <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = PUDP, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("PUDP") +
  scale_fill_brewer(palette = "Dark2")

KDM6A_ad <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = KDM6A, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("KDM6A") +
  scale_fill_brewer(palette = "Dark2")

SYAP1_ad <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = SYAP1, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("SYAP1") +
  scale_fill_brewer(palette = "Dark2")

JPX_ad <- ggplot(pheno_snc_rosmap.data_expr_cell.data, aes(x = cell, y = JPX, fill = FK_final_CaseControl)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T, scale = "width") +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE, outlier.shape = NA) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  guides(fill=guide_legend(title="AD status")) +
  ylab("Log Normalized Expression") +
  ggtitle("JPX") +
  scale_fill_brewer(palette = "Dark2")


plots <- list(KDM5C,KDM5C_ad,PUDP,PUDP_ad,KDM6A,KDM6A_ad,SYAP1,SYAP1_ad,JPX,JPX_ad)
all_genes_single <- plot_grid(plotlist = plots, ncol = 2)
ggsave(all_genes_single, file = "figures/all_genes_single.png", width = 5000, height = 5000, units = "px")


### SEM
mediation_rosmap.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov8 + InferredCov10 + InferredCov18 + InferredCov20 + InferredCov24 + InferredCov25 + InferredCov27 + InferredCov29 + InferredCov34 + InferredCov35 + InferredCov36 + InferredCov41 + InferredCov47 + InferredCov49 + InferredCov51 + InferredCov56 + InferredCov60
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov8 + InferredCov10 + InferredCov18 + InferredCov20 + InferredCov24 + InferredCov25 + InferredCov27 + InferredCov29 + InferredCov34 + InferredCov35 + InferredCov36 + InferredCov41 + InferredCov47 + InferredCov49 + InferredCov51 + InferredCov56 + InferredCov60
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

rosmap_transcripts <- c("KDM5C","PUDP","KDM6A","SYAP1","JPX")

# Astrocytes
pheno_snc_rosmap_ast.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Ast", ]

mediation_rosmap_ast.list <- mclapply(pheno_snc_rosmap_ast.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_ast.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_ast.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_ast.data <- do.call(rbind, mediation_rosmap_ast.list)
mediation_rosmap_ast.data$cell <- "Astrocytes"

# Cell adhesion molecules
pheno_snc_rosmap_cam.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "CAMs", ]

mediation_rosmap_cam.list <- mclapply(pheno_snc_rosmap_cam.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_cam.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_cam.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_cam.data <- do.call(rbind, mediation_rosmap_cam.list)
mediation_rosmap_cam.data$cell <- "Cell Adhesion Molecules"

# Endothelial
pheno_snc_rosmap_end.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "End", ]

mediation_rosmap_end.list <- mclapply(pheno_snc_rosmap_end.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_end.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_end.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_end.data <- do.call(rbind, mediation_rosmap_end.list)
mediation_rosmap_end.data$cell <- "Endothelial"

# Excitatory neurons
pheno_snc_rosmap_exc.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Exc", ]

mediation_rosmap_exc.list <- mclapply(pheno_snc_rosmap_exc.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_exc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_exc.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_exc.data <- do.call(rbind, mediation_rosmap_exc.list)
mediation_rosmap_exc.data$cell <- "Excitatory Neurons"

# Fibroblast
pheno_snc_rosmap_fib.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Fib", ]

mediation_rosmap_fib.list <- mclapply(pheno_snc_rosmap_fib.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_fib.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_fib.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_fib.data <- do.call(rbind, mediation_rosmap_fib.list)
mediation_rosmap_fib.data$cell <- "Fibroblast"

# Inh
pheno_snc_rosmap_inh.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Inh", ]

mediation_rosmap_inh.list <- mclapply(pheno_snc_rosmap_inh.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_inh.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_inh.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_inh.data <- do.call(rbind, mediation_rosmap_inh.list)
mediation_rosmap_inh.data$cell <- "Inhibitory Neurons"

# Mic
pheno_snc_rosmap_mic.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Mic", ]

mediation_rosmap_mic.list <- mclapply(pheno_snc_rosmap_mic.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_mic.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_mic.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_mic.data <- do.call(rbind, mediation_rosmap_mic.list)
mediation_rosmap_mic.data$cell <- "Microglia"

# Oligodendrocyte
pheno_snc_rosmap_oli.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Oli", ]

mediation_rosmap_oli.list <- mclapply(pheno_snc_rosmap_oli.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_oli.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_oli.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_oli.data <- do.call(rbind, mediation_rosmap_oli.list)
mediation_rosmap_oli.data$cell <- "Oligodendrocytes"

# Oligodendrocyte precursor cells 
pheno_snc_rosmap_opc.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "OPC", ]

mediation_rosmap_opc.list <- mclapply(pheno_snc_rosmap_opc.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_opc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_opc.data, estimator = "WLSMV")
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
}, mc.cores = 8)
mediation_rosmap_opc.data <- do.call(rbind, mediation_rosmap_opc.list)
mediation_rosmap_opc.data$cell <- "Oligodendrocyte Preursor Cells"

# Pericyte
pheno_snc_rosmap_per.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "Per", ]

mediation_rosmap_per.list <- mclapply(pheno_snc_rosmap_per.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_per.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_per.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_per.data <- do.call(rbind, mediation_rosmap_per.list)
mediation_rosmap_per.data$cell <- "Pericytes"

# Smooth Muscle Cells
pheno_snc_rosmap_smc.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "SMC", ]

mediation_rosmap_smc.list <- mclapply(pheno_snc_rosmap_smc.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_smc.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_smc.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_smc.data <- do.call(rbind, mediation_rosmap_smc.list)
mediation_rosmap_smc.data$cell <- "Smooth Muscle Cells"

# T Cells
pheno_snc_rosmap_T.data <- pheno_snc_rosmap.data[pheno_snc_rosmap.data$cell == "T", ]

mediation_rosmap_T.list <- mclapply(pheno_snc_rosmap_T.data[,..rosmap_transcripts], function(transcript){
  pheno_snc_rosmap_T.data$transcript <- scale(transcript)
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_snc_rosmap_T.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
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
}, mc.cores = 8)
mediation_rosmap_T.data <- do.call(rbind, mediation_rosmap_T.list)
mediation_rosmap_T.data$cell <- "T Cells"

mediation_rosmap_snRNAseq.data <- rbind(mediation_rosmap_ast.data,mediation_rosmap_exc.data,mediation_rosmap_inh.data,mediation_rosmap_mic.data,mediation_rosmap_oli.data,mediation_rosmap_opc.data)
write.csv(mediation_rosmap_snRNAseq.data, file = "results/mediation_rosmap_snRNAseq.csv")
