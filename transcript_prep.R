library(data.table)

###### Preparation of expression matrices for further analysis

### Autosomes and X-chromosome
## Read in expression matrix, split by chromosome and transpose, so that rows participants and columns genes
# ROSMAP
rosmap_dlpfc_chr <- fread("expression_matrices_and_PEER/ROSMAP_DLPFC.autosomesX.allSamples_n636.14April2025.expression.bed.gz", data.table =  F)
rosmap_dlpfc.list <- split(rosmap_dlpfc_chr, rosmap_dlpfc_chr$'#chr')

lapply(rosmap_dlpfc.list, function(rosmap_dlpfc) {
  row.names(rosmap_dlpfc) <- rosmap_dlpfc$gene_id
  filename <- paste0(paste("transcriptome/ROSMAP/rosmap_dlpfc",rosmap_dlpfc$`#chr`[1],sep="_"),".txt")
  rosmap_dlpfc <- rosmap_dlpfc[,5:dim(rosmap_dlpfc)[2]]
  rosmap_dlpfc_t <- t(rosmap_dlpfc)
  rosmap_dlpfc.mat <- cbind(row.names(rosmap_dlpfc_t), rosmap_dlpfc_t)
  rosmap_dlpfc.data <- as.data.table(rosmap_dlpfc.mat)
  setnames(rosmap_dlpfc.data,"V1","IID")
  fwrite(rosmap_dlpfc.data, file = filename, sep = " ")
})

# Mayo
mayo_tcx_chr <- fread("expression_matrices_and_PEER/MayoRNASeq_TCX.autosomesX.allSamples_n264.14April2025.expression.bed.gz", data.table =  F)
mayo_tcx.list <- split(mayo_tcx_chr, mayo_tcx_chr$'#chr')

lapply(mayo_tcx.list, function(mayo_tcx) {
  row.names(mayo_tcx) <- mayo_tcx$gene_id
  filename <- paste0(paste("transcriptome/Mayo/mayo_tcx",mayo_tcx$`#chr`[1],sep="_"),".txt")
  mayo_tcx <- mayo_tcx[,5:dim(mayo_tcx)[2]]
  mayo_tcx_t <- t(mayo_tcx)
  mayo_tcx.mat <- cbind(row.names(mayo_tcx_t), mayo_tcx_t)
  mayo_tcx.data <- as.data.table(mayo_tcx.mat)
  setnames(mayo_tcx.data,"V1","IID")
  fwrite(mayo_tcx.data, file = filename, sep = " ")
})

# MSBB BA10
msbb_ba10_chr <- fread("expression_matrices_and_PEER/MSBB_BA10.autosomesX.allSamples_n227.14April2025.expression.bed.gz", data.table =  F)
msbb_ba10.list <- split(msbb_ba10_chr, msbb_ba10_chr$'#chr')

lapply(msbb_ba10.list, function(msbb_ba10) {
  row.names(msbb_ba10) <- msbb_ba10$gene_id
  filename <- paste0(paste("transcriptome/MSBB_BA10/msbb_ba10",msbb_ba10$`#chr`[1],sep="_"),".txt")
  msbb_ba10 <- msbb_ba10[,5:dim(msbb_ba10)[2]]
  msbb_ba10_t <- t(msbb_ba10)
  msbb_ba10.mat <- cbind(row.names(msbb_ba10_t), msbb_ba10_t)
  msbb_ba10.data <- as.data.table(msbb_ba10.mat)
  setnames(msbb_ba10.data,"V1","IID")
  fwrite(msbb_ba10.data, file = filename, sep = " ")
})

# MSBB BA22
msbb_ba22_chr <- fread("expression_matrices_and_PEER/MSBB_BA22.autosomesX.allSamples_n202.14April2025.expression.bed.gz", data.table =  F)
msbb_ba22.list <- split(msbb_ba22_chr, msbb_ba22_chr$'#chr')

lapply(msbb_ba22.list, function(msbb_ba22) {
  row.names(msbb_ba22) <- msbb_ba22$gene_id
  filename <- paste0(paste("transcriptome/MSBB_BA22/msbb_ba22",msbb_ba22$`#chr`[1],sep="_"),".txt")
  msbb_ba22 <- msbb_ba22[,5:dim(msbb_ba22)[2]]
  msbb_ba22_t <- t(msbb_ba22)
  msbb_ba22.mat <- cbind(row.names(msbb_ba22_t), msbb_ba22_t)
  msbb_ba22.data <- as.data.table(msbb_ba22.mat)
  setnames(msbb_ba22.data,"V1","IID")
  fwrite(msbb_ba22.data, file = filename, sep = " ")
})

# MSBB BA36
msbb_ba36_chr <- fread("expression_matrices_and_PEER/MSBB_BA36.autosomesX.allSamples_n180.14April2025.expression.bed.gz", data.table =  F)
msbb_ba36.list <- split(msbb_ba36_chr, msbb_ba36_chr$'#chr')

lapply(msbb_ba36.list, function(msbb_ba36) {
  row.names(msbb_ba36) <- msbb_ba36$gene_id
  filename <- paste0(paste("transcriptome/MSBB_BA36/msbb_ba36",msbb_ba36$`#chr`[1],sep="_"),".txt")
  msbb_ba36 <- msbb_ba36[,5:dim(msbb_ba36)[2]]
  msbb_ba36_t <- t(msbb_ba36)
  msbb_ba36.mat <- cbind(row.names(msbb_ba36_t), msbb_ba36_t)
  msbb_ba36.data <- as.data.table(msbb_ba36.mat)
  setnames(msbb_ba36.data,"V1","IID")
  fwrite(msbb_ba36.data, file = filename, sep = " ")
})

# MSBB BA44
msbb_ba44_chr <- fread("expression_matrices_and_PEER/MSBB_BA44.autosomesX.allSamples_n193.14April2025.expression.bed.gz", data.table =  F)
msbb_ba44.list <- split(msbb_ba44_chr, msbb_ba44_chr$'#chr')

lapply(msbb_ba44.list, function(msbb_ba44) {
  row.names(msbb_ba44) <- msbb_ba44$gene_id
  filename <- paste0(paste("transcriptome/MSBB_BA44/msbb_ba44",msbb_ba44$`#chr`[1],sep="_"),".txt")
  msbb_ba44 <- msbb_ba44[,5:dim(msbb_ba44)[2]]
  msbb_ba44_t <- t(msbb_ba44)
  msbb_ba44.mat <- cbind(row.names(msbb_ba44_t), msbb_ba44_t)
  msbb_ba44.data <- as.data.table(msbb_ba44.mat)
  setnames(msbb_ba44.data,"V1","IID")
  fwrite(msbb_ba44.data, file = filename, sep = " ")
})

# Merged ROSMAP, Mayo, MSBB BA10
amp_ad_chr <- fread("expression_matrices_and_PEER/autosomesX.rosmapDLPFC_mayoTCX_msbbBA10_innerJoint_normalizedSeparately_allSamples_14April2025.bed.gz", data.table =  F)
amp_ad.list <- split(amp_ad_chr, amp_ad_chr$'#chr')

lapply(amp_ad.list, function(amp_ad) {
  row.names(amp_ad) <- amp_ad$gene_id
  filename <- paste0(paste("transcriptome/amp_ad/amp_ad",amp_ad$`#chr`[1],sep="_"),".txt")
  amp_ad <- amp_ad[,5:dim(amp_ad)[2]]
  amp_ad_t <- t(amp_ad)
  amp_ad.mat <- cbind(row.names(amp_ad_t), amp_ad_t)
  amp_ad.data <- as.data.table(amp_ad.mat)
  setnames(amp_ad.data,"V1","IID")
  fwrite(amp_ad.data, file = filename, sep = " ")
})

# Y-chromosome
# ROSMAP
rosmap_dlpfc_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.ROSMAP_DLPFC.autosomesXY.allSamples_n229.14April2025.expression.bed.gz", data.table =  F)
rosmap_dlpfc.list <- split(rosmap_dlpfc_chr, rosmap_dlpfc_chr$'#chr')

rosmap_dlpfc <- rosmap_dlpfc.list[[24]]
row.names(rosmap_dlpfc) <- rosmap_dlpfc$gene_id
filename <- paste0(paste("transcriptome/ROSMAP/rosmap_dlpfc",rosmap_dlpfc$`#chr`[24],sep="_"),".txt")
rosmap_dlpfc <- rosmap_dlpfc[,5:dim(rosmap_dlpfc)[2]]
rosmap_dlpfc_t <- t(rosmap_dlpfc)
rosmap_dlpfc.mat <- cbind(row.names(rosmap_dlpfc_t), rosmap_dlpfc_t)
rosmap_dlpfc.data <- as.data.table(rosmap_dlpfc.mat)
setnames(rosmap_dlpfc.data,"V1","IID")
fwrite(rosmap_dlpfc.data, file = filename, sep = " ")

# Mayo
mayo_tcx_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.MayoRNASeq_TCX.autosomesXY.allSamples_n129.14April2025.expression.bed.gz", data.table =  F)
mayo_tcx.list <- split(mayo_tcx_chr, mayo_tcx_chr$'#chr')

mayo_tcx <- mayo_tcx.list[[24]]
row.names(mayo_tcx) <- mayo_tcx$gene_id
filename <- paste0(paste("transcriptome/Mayo/mayo_tcx",mayo_tcx$`#chr`[24],sep="_"),".txt")
mayo_tcx <- mayo_tcx[,5:dim(mayo_tcx)[2]]
mayo_tcx_t <- t(mayo_tcx)
mayo_tcx.mat <- cbind(row.names(mayo_tcx_t), mayo_tcx_t)
mayo_tcx.data <- as.data.table(mayo_tcx.mat)
setnames(mayo_tcx.data,"V1","IID")
fwrite(mayo_tcx.data, file = filename, sep = " ")

# MSBB BA10
msbb_ba10_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.MSBB_BA10.autosomesXY.allSamples_n81.14April2025.expression.bed.gz", data.table =  F)
msbb_ba10.list <- split(msbb_ba10_chr, msbb_ba10_chr$'#chr')

msbb_ba10 <- msbb_ba10.list[[24]]
row.names(msbb_ba10) <- msbb_ba10$gene_id
filename <- paste0(paste("transcriptome/MSBB_BA10/msbb_ba10",msbb_ba10$`#chr`[24],sep="_"),".txt")
msbb_ba10 <- msbb_ba10[,5:dim(msbb_ba10)[2]]
msbb_ba10_t <- t(msbb_ba10)
msbb_ba10.mat <- cbind(row.names(msbb_ba10_t), msbb_ba10_t)
msbb_ba10.data <- as.data.table(msbb_ba10.mat)
setnames(msbb_ba10.data,"V1","IID")
fwrite(msbb_ba10.data, file = filename, sep = " ")

# MSBB BA22
msbb_ba22_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.MSBB_BA22.autosomesXY.allSamples_n70.14April2025.expression.bed.gz", data.table =  F)
msbb_ba22.list <- split(msbb_ba22_chr, msbb_ba22_chr$'#chr')

msbb_ba22 <- msbb_ba22.list[[24]]
row.names(msbb_ba22) <- msbb_ba22$gene_id
filename <- paste0(paste("transcriptome/MSBB_BA22/msbb_ba22",msbb_ba22$`#chr`[24],sep="_"),".txt")
msbb_ba22 <- msbb_ba22[,5:dim(msbb_ba22)[2]]
msbb_ba22_t <- t(msbb_ba22)
msbb_ba22.mat <- cbind(row.names(msbb_ba22_t), msbb_ba22_t)
msbb_ba22.data <- as.data.table(msbb_ba22.mat)
setnames(msbb_ba22.data,"V1","IID")
fwrite(msbb_ba22.data, file = filename, sep = " ")

# MSBB BA36
msbb_ba36_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.MSBB_BA36.autosomesXY.allSamples_n68.14April2025.expression.bed.gz", data.table =  F)
msbb_ba36.list <- split(msbb_ba36_chr, msbb_ba36_chr$'#chr')

msbb_ba36 <- msbb_ba36.list[[24]]
row.names(msbb_ba36) <- msbb_ba36$gene_id
filename <- paste0(paste("transcriptome/MSBB_BA36/msbb_ba36",msbb_ba36$`#chr`[24],sep="_"),".txt")
msbb_ba36 <- msbb_ba36[,5:dim(msbb_ba36)[2]]
msbb_ba36_t <- t(msbb_ba36)
msbb_ba36.mat <- cbind(row.names(msbb_ba36_t), msbb_ba36_t)
msbb_ba36.data <- as.data.table(msbb_ba36.mat)
setnames(msbb_ba36.data,"V1","IID")
fwrite(msbb_ba36.data, file = filename, sep = " ")

# MSBB BA44
msbb_ba44_chr <- fread("expression_matrices_and_PEER/maleOnly/malesOnly.MSBB_BA44.autosomesXY.allSamples_n68.14April2025.expression.bed.gz", data.table =  F)
msbb_ba44.list <- split(msbb_ba44_chr, msbb_ba44_chr$'#chr')

msbb_ba44 <- msbb_ba44.list[[24]]
row.names(msbb_ba44) <- msbb_ba44$gene_id
filename <- paste0(paste("transcriptome/MSBB_BA44/msbb_ba44",msbb_ba44$`#chr`[24],sep="_"),".txt")
msbb_ba44 <- msbb_ba44[,5:dim(msbb_ba44)[2]]
msbb_ba44_t <- t(msbb_ba44)
msbb_ba44.mat <- cbind(row.names(msbb_ba44_t), msbb_ba44_t)
msbb_ba44.data <- as.data.table(msbb_ba44.mat)
setnames(msbb_ba44.data,"V1","IID")
fwrite(msbb_ba44.data, file = filename, sep = " ")

