library(pacman)
p_load(data.table, bsseq, ggplot2, DSS, data.table, annotatr, GenomicRanges, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, pheatmap)
options(scipen=999)
outdatadir <- "./outputs/"

# MAIN SCRIPT
sample_info <- readRDS("outputs/sample_info.RDS")  
setDT(sample_info)

if(!file.exists(paste0(outdatadir, "bsseq_object.RDS"))) {
  methylation_dt <- readRDS("outputs/methylation_dt.RDS")
  coverage_dt <- readRDS("outputs/coverage_dt.RDS")
  
  sample_cols <- names(methylation_dt)[!names(methylation_dt) %in% c("cpg_id", "chr", "position")]

  # Create BSseq object
  chr <- methylation_dt$chr
  pos <- methylation_dt$position

  M_matrix <- matrix(NA, nrow = nrow(methylation_dt), ncol = length(sample_cols))
  Cov_matrix <- matrix(NA, nrow = nrow(methylation_dt), ncol = length(sample_cols))

  for(i in 1:length(sample_cols)) {
    sample <- sample_cols[i]
    
    meth <- methylation_dt[[sample]]
    cov <- coverage_dt[[sample]]
    
    M_matrix[, i] <- round(meth * cov, 0)
    Cov_matrix[, i] <- cov
    
    M_matrix[is.na(M_matrix[, i]), i] <- 0
    Cov_matrix[is.na(Cov_matrix[, i]), i] <- 0
  }

  colnames(M_matrix) <- sample_cols
  colnames(Cov_matrix) <- sample_cols

  BS <- BSseq(chr = chr,
              pos = pos,
              M = M_matrix,
              Cov = Cov_matrix,
              sampleNames = sample_cols)

  pData(BS) <- sample_info
  rownames(pData(BS)) <- sample_info$sample_id
  saveRDS(BS, paste0(outdatadir, "bsseq_object.RDS")) } else {
  BS <- readRDS(paste0(outdatadir, "bsseq_object.RDS"))
  } 

# Define comparison
group1_samples <- sample_info[tissue == "T1" & condition == "healthy", sample_id]
group2_samples <- sample_info[tissue == "T1" & condition == "tumor", sample_id]

# Filtering BS data. What can happen if we skip this? 

cov_matrix <- getCoverage(BS)
samples_covered <- rowSums(cov_matrix[, 1:4] >= 10, na.rm = TRUE)
keep_cpgs <- samples_covered == 4
BS_filtered <- BS[keep_cpgs]
saveRDS(BS_filtered, paste0(outdatadir, "BS_filtered.RDS"))

# run the test
if(!(file.exists(paste0(outdatadir, "dss_dmlTest_result.RDS")))) {
  dmlTest <- DMLtest(BS_filtered, 
                    group1 = group1_samples,
                    group2 = group2_samples,
                    ncores = 8, 
                    smoothing = FALSE)
  saveRDS(dmlTest, paste0(outdatadir, "dss_dmlTest_result.RDS"))} else {
    dmlTest <- readRDS(paste0(outdatadir, "dss_dmlTest_result.RDS"))
}

# Call DMCs p < 0.05, delta > 0.1)

dmcs_delta <- callDML(dmlTest, 
                p.threshold = 0.05,
                delta = 0.1)
# Summary
dmcs_delta_dt <- as.data.table(dmcs_delta)
n_hyper_delta <- sum(dmcs_delta_dt$diff > 0)
n_hypo_delta <- sum(dmcs_delta_dt$diff < 0)


dmcs_fdr <- dmlTest[dmlTest$fdr < 0.05, ]
dmcs_fdr_dt <- as.data.table(dmcs_fdr)
n_hyper_fdr <- sum(dmcs_fdr_dt$diff > 0)
n_hypo_fdr <- sum(dmcs_fdr_dt$diff < 0)

dim(dmlTest)
dim(dmcs_delta_dt)
dim(dmcs_fdr_dt)

# https://chatgpt.com/c/6a072a96-7a48-8332-9389-38c3552a3389

# =========================================================
# 4. DSS vs Wilcoxon comparison
# =========================================================

# ---------------------------------------------------------
# Wilcoxon test on chr1
# ---------------------------------------------------------

meth <- getMeth(BS_filtered, type = "raw")

wilcox_p <- apply(meth, 1, function(x) {

  wilcox.test(
    x[group1_samples],
    x[group2_samples]
  )$p.value

})

wilcox_fdr <- p.adjust(
  wilcox_p,
  method = "fdr"
)

group1_mean <- rowMeans(
  meth[, group1_samples],
  na.rm = TRUE
)

group2_mean <- rowMeans(
  meth[, group2_samples],
  na.rm = TRUE
)

delta <- group1_mean - group2_mean

wilcox_dmc_analysis <- data.table(
  chr = as.character(seqnames(BS_filtered)),
  pos = start(BS_filtered),
  fdr = wilcox_fdr,
  p_value = wilcox_p,
  diff = delta
)

wilcox_dmcs <- wilcox_dmc_analysis[p_value < 0.5 & abs(diff) > 0.1]

# =========================================================
# Create Venn diagram
# =========================================================
library(VennDiagram)

beta_binomial_sites <- paste(
  dmcs_delta_dt$chr,
  dmcs_delta_dt$pos
)

wilcox_sites <- paste(
  wilcox_dmcs$chr,
  wilcox_dmcs$pos
)

venn.diagram(
  x = list(
    Beta_binomial = beta_binomial_sites,
    Wilcoxon = wilcox_sites
  ),
  filename = paste0(
    outdatadir,
    "BetaBin_vs_Wilcoxon_venn.tif"
  )
)

# =========================================================
# Compare methods
# =========================================================

comparison <- data.table(
  Method = c(
    "Beta_binomial",
    "Wilcoxon"
  ),

  Num_DMCs = c(
    nrow(dmcs_delta_dt),
    nrow(wilcox_dmcs)
  ),

  Mean_Effect_Size = c(
    mean(abs(dmcs_delta_dt$diff), na.rm = TRUE),
    mean(abs(wilcox_dmcs$diff), na.rm = TRUE)
  )
)

print(comparison)

write.csv(
  comparison,
  paste0(outdatadir, "DSS_vs_Wilcoxon_summary.csv"),
  row.names = FALSE
)

# =========================================================
# 5. Threshold exploration
# =========================================================

res1 <- callDML(
  dmlTest,
  p.threshold = 0.05,
  delta = 0.1
)

res2 <- callDML(
  dmlTest,
  p.threshold = 0.01,
  delta = 0.1
)

res3 <- callDML(
  dmlTest,
  p.threshold = 0.05,
  delta = 0.2
)

res4 <- callDML(
  dmlTest,
  p.threshold = 0.01,
  delta = 0.2
)

# =========================================================
# Summarize threshold exploration
# =========================================================

threshold_summary <- data.table(

  Threshold = c(
    "p<0.05, delta>0.1",
    "p<0.01, delta>0.1",
    "p<0.05, delta>0.2",
    "p<0.01, delta>0.2"
  ),

  Num_DMCs = c(
    nrow(res1),
    nrow(res2),
    nrow(res3),
    nrow(res4)
  )
)

print(threshold_summary)

write.csv(
  threshold_summary,
  paste0(outdatadir, "threshold_exploration_summary.csv"),
  row.names = FALSE
)
