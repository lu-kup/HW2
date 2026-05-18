library(pacman)
p_load(data.table, bsseq, ggplot2, DSS, data.table, annotatr, GenomicRanges, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, pheatmap)
options(scipen=999)
outdatadir <- "./outputs/"

# Identification script
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
n_hypo_delta <- sum(dmcs_delta_dt$diff > 0)
n_hyper_delta <- sum(dmcs_delta_dt$diff < 0)


dmcs_fdr <- dmlTest[dmlTest$fdr < 0.05, ]
dmcs_fdr_dt <- as.data.table(dmcs_fdr)
n_hypo_fdr <- sum(dmcs_fdr_dt$diff > 0)
n_hyper_fdr <- sum(dmcs_fdr_dt$diff < 0)

dim(dmlTest)
dim(dmcs_delta_dt)
dim(dmcs_fdr_dt)

# =========================================================
# 4. DSS vs Wilcoxon comparison
# =========================================================

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


# Annotation

if(!file.exists(paste0(outdatadir, "annotation_hg38.RDS"))) {
  annots <- c('hg38_genes_promoters')
  annotations <- build_annotations(genome = 'hg38', annotations = annots)
  saveRDS(annotations, paste0(outdatadir, "annotation_hg38.RDS")) } else {
  annotations <- readRDS(paste0(outdatadir, "annotation_hg38.RDS"))
  }

dmcs_gr <- GRanges(
  seqnames = dmcs_delta_dt$chr,
  ranges = IRanges(start = dmcs_delta_dt$pos, width = 1),
  diff = dmcs_delta_dt$diff,
  pval = dmcs_delta_dt$pval,
  fdr = dmcs_delta_dt$fdr,
  stat = dmcs_delta_dt$stat
)

dmcs_annotated <- annotate_regions(
  regions = dmcs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)

dmcs_annot_dt <- as.data.table(dmcs_annotated)
dmcs_annot_dt[, annot_type := gsub("hg38_", "", annot.type)]
dmcs_annot_dt[, annot_type := gsub("genes_", "", annot_type)]
dmcs_annot_dt[, annot_type := gsub("cpg_", "", annot_type)]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)
nearest_genes <- distanceToNearest(dmcs_gr, all_genes)

dmcs_nearest <- data.table(
  dmc_index = queryHits(nearest_genes),
  gene_id = names(all_genes)[subjectHits(nearest_genes)],
  distance_to_gene = mcols(nearest_genes)$distance
)

gene_symbols <- select(org.Hs.eg.db,
                      keys = unique(dmcs_nearest$gene_id),
                      columns = c("SYMBOL", "GENENAME"),
                      keytype = "ENTREZID")

gene_symbols <- as.data.table(gene_symbols)
gene_symbols <- gene_symbols[!duplicated(ENTREZID)]

dmcs_nearest <- merge(dmcs_nearest, gene_symbols, 
                     by.x = "gene_id", by.y = "ENTREZID",
                     all.x = TRUE)

setnames(dmcs_nearest, "SYMBOL", "nearest_gene_symbol")
setnames(dmcs_nearest, "GENENAME", "nearest_gene_name")

dmcs_with_genes <- copy(dmcs_delta_dt)
dmcs_with_genes[, dmc_index := 1:.N]

dmcs_with_genes <- merge(dmcs_with_genes, dmcs_nearest,
                         by = "dmc_index", all.x = TRUE)

# Take the first annotation for each DMC
dmcs_features <- dmcs_annot_dt[, .SD[1], by = .(seqnames, start)]

setnames(dmcs_features, c("seqnames", "start"), c("chr", "pos"))

dmcs_final <- merge(dmcs_with_genes, 
                   dmcs_features[, .(chr, pos, annot_type)],
                   by = c("chr", "pos"), all.x = TRUE)

dmcs_final[, direction := ifelse(diff < 0, "Hypermethylated", "Hypomethylated")]

# Reorder columns
setcolorder(dmcs_final, c("chr", "pos", "diff", "pval", "fdr", "stat",
                         "direction", "annot_type", 
                         "gene_id", "nearest_gene_symbol", "nearest_gene_name",
                         "distance_to_gene"))

saveRDS(dmcs_annotated, paste0(outdatadir, "dmcs_annotated_granges.RDS"))
saveRDS(dmcs_final, paste0(outdatadir, "dmcs_final.RDS"))

promoter_dmcs <- dmcs_final[grepl("promoter", annot_type, ignore.case = TRUE)]

write.csv(
  summary(promoter_dmcs),
  file = "outputs/promoter_dmcs_summary.csv",
  row.names = FALSE
)

write.csv(
  summary(dmcs_final),
  file = "outputs/dmcs_final_summary.csv",
  row.names = FALSE
)

# Enrichment

feature_counts <- dmcs_final[, .N, by = annot_type]
feature_counts[, percentage := 100 * N / sum(N)]

dmcs_final[, feature_category := fcase(
  grepl("promoter", annot_type, ignore.case = TRUE), "Promoters",
  default = "Other"
)]

feature_summary <- dmcs_final[, .N, by = feature_category]
feature_summary[, percentage := 100 * N / sum(N)]

all_cpgs <- BS_filtered
all_cpgs_gr <- granges(all_cpgs)

background_annotated <- annotate_regions(
  regions = all_cpgs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE
)

background_dt <- as.data.table(background_annotated)
background_dt[, annot_type := gsub("hg19_", "", annot.type)]
background_dt[, annot_type := gsub("genes_", "", annot_type)]

background_features <- background_dt[, .SD[1], by = .(seqnames, start)]
setnames(background_features, c("seqnames", "start"), c("chr", "pos"))

all_cpgs_dt <- data.table(
  chr = as.character(seqnames(BS_filtered)),
  pos = start(BS_filtered)
)

background_final <- merge(all_cpgs_dt, 
                   background_features[, .(chr, pos, annot_type)],
                   by = c("chr", "pos"), all.x = TRUE)

background_final[, feature_category := fcase(
  grepl("promoter", annot_type, ignore.case = TRUE), "Promoters",
  default = "Other"
)]

background_summary <- background_final[, .N, by = feature_category]
background_summary[, bg_percentage := 100 * N / sum(N)]

# Merge and calculate enrichment
enrichment_table <- merge(feature_summary, 
                         background_summary[, .(feature_category, bg_percentage)],
                         by = "feature_category")

enrichment_table[, enrichment := percentage / bg_percentage]
enrichment_table[, log2_enrichment := log2(enrichment)]

# Diagrams

ggplot(feature_summary, aes(x = reorder(feature_category, -percentage), 
                                   y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.5, size = 3.5) +
  labs(title = "DMC Distribution across Genomic Features",
       x = "Genomic Feature",
       y = "Percentage of DMCs (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "diagrams/dmc_genomic_features.png",
  width = 8,
  height = 6,
  dpi = 300
)


ggplot(enrichment_table, aes(x = reorder(feature_category, enrichment), 
                                    y = enrichment)) +
  geom_bar(stat = "identity", aes(fill = enrichment > 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                    labels = c("Depleted", "Enriched")) +
  geom_text(aes(label = round(enrichment, 2)), hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "DMC Enrichment in Genomic Features",
       x = "Genomic Feature",
       y = "Enrichment (Observed / Expected)",
       fill = "") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  filename = "diagrams/dmc_enrichment.png",
  width = 8,
  height = 6,
  dpi = 300
)


# Features by direction

feature_by_direction <- dmcs_final[, .N, by = .(feature_category, direction)]

feature_by_direction[, total := sum(N), by = direction]
feature_by_direction[, percentage := 100 * N / total]


ggplot(feature_by_direction, aes(x = feature_category, y = percentage, 
                                        fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Hypermethylated" = "#E64B35",
                               "Hypomethylated" = "#4DBBD5")) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = "Genomic Feature Preferences: Hyper vs Hypo DMCs",
       x = "Genomic Feature",
       y = "Percentage of DMCs (%)",
       fill = "Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(
  filename = "diagrams/dmc_feature_by_direction.png",
  width = 8,
  height = 6,
  dpi = 300
)


# Volcano plot

volcano_data <- copy(dmcs_final)
volcano_data[, neg_log_fdr := -log10(fdr)]

p <- ggplot(volcano_data, aes(x = diff, y = neg_log_fdr, color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Hypermethylated" = "#E64B35",
                                "Hypomethylated" = "#4DBBD5")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot - DMCs",
       x = "Methylation Difference (Δβ)",
       y = "-log10(FDR)",
       color = "Direction") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  filename = "outputs/volcano_plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# MA-equivalent plot

ma_data <- copy(dmcs_final)
ma_data[, avg_meth := (mu1 + mu2) / 2]

p <- ggplot(ma_data, aes(x = avg_meth, y = diff, color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Hypermethylated" = "#E64B35",
                                "Hypomethylated" = "#4DBBD5")) +
  labs(title = "MA-equivalent plot - DMCs",
       x = "Average Methylation",
       y = "Methylation Difference (Δβ)",
       color = "Direction") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  filename = "outputs/ma_plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# DMR calling

dmrs <- callDMR(dmlTest, minCG=3, minlen=50, p.threshold=0.05)
dmrs2 <- callDMR(dmlTest, minCG=5, minlen=100, p.threshold=0.01)
dmrs3 <- callDMR(dmlTest, minCG=7, minlen=200, p.threshold=0.01)

write.csv(
  summary(dmrs),
  file = "outputs/dmrs1_summary.csv",
  row.names = FALSE
)

write.csv(
  summary(dmrs2),
  file = "outputs/dmrs2_summary.csv",
  row.names = FALSE
)

write.csv(
  summary(dmrs3),
  file = "outputs/dmrs3_summary.csv",
  row.names = FALSE
)

# Characterize your DMRs

library(dplyr)

dmrs_char_strict <- dmrs2 %>%
  mutate(
    diff = meanMethy1 - meanMethy2,
    direction = ifelse(diff < 0,
                       "Hypermethylated",
                       "Hypomethylated"),
    abs_diff = abs(diff)
  )

dmrs_char <- dmrs %>%
  mutate(
    diff = meanMethy1 - meanMethy2,
    direction = ifelse(diff < 0,
                       "Hypermethylated",
                       "Hypomethylated"),
    abs_diff = abs(diff)
  )


p <- ggplot(dmrs_char, aes(x = length)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(breaks = seq(0, 10000, by = 1000)) +
  labs(
    title = "DMR Length Distribution",
    x = "DMR Length (bp)",
    y = "Count"
  ) +
  theme_bw()

ggsave("diagrams/dmr_length_distribution.png", width = 7, height = 5, dpi = 300)


ggplot(dmrs_char, aes(x = nCG)) +
  geom_histogram(binwidth = 1, boundary = 0.5) +
  scale_x_continuous(
    breaks = seq(0, max(dmrs$nCG, na.rm = TRUE), by = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 1000, by = 10)) +
  labs(
    title = "CpGs per DMR Distribution",
    x = "Number of CpGs per DMR",
    y = "Count"
  ) +
  theme_bw()

ggsave("diagrams/dmr_cpg_distribution.png", width = 7, height = 5, dpi = 300)


ggplot(dmrs_char_strict, aes(x = direction, fill = direction)) +
  geom_bar() +
  labs(
    title = "Hypermethylated vs Hypo (minCG = 5, minlen = 100, p = 0.01)",
    x = "",
    y = "Count"
  ) +
  theme_bw()

ggsave("diagrams/hyper_vs_hypo_dmrs_strict.png",
       width = 6,
       height = 5,
       dpi = 300)


ggplot(dmrs_char, aes(x = direction, fill = direction)) +
  geom_bar() +
  labs(
    title = "Hypermethylated vs Hypo (minCG = 3, minlen = 50, p = 0.05)",
    x = "",
    y = "Count"
  ) +
  theme_bw()

ggsave("diagrams/hyper_vs_hypo_dmrs.png",
       width = 6,
       height = 5,
       dpi = 300)


ggplot(dmrs_char, aes(x = start, color = chr)) +
  geom_density() +
  labs(
    title = "DMR Density Across Genome",
    x = "Genomic Position",
    y = "Density"
  ) +
  theme_bw()

ggsave("diagrams/dmr_density_genome.png",
       width = 8,
       height = 5,
       dpi = 300)


ggplot(dmrs_char, aes(x = length, y = abs_diff)) +
  scale_x_continuous(breaks = seq(0, 10000, by = 1000)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "DMR Length vs Effect Size",
    x = "DMR Length (bp)",
    y = "Absolute Methylation Difference"
  ) +
  theme_bw()

ggsave("diagrams/dmr_length_vs_effect.png",
       width = 7,
       height = 5,
       dpi = 300)

saveRDS(dmrs_char, paste0(outdatadir, "dmrs_char.RDS"))

# Compare DMCs and DMRs

library(GenomicRanges)

dmc_gr <- GRanges(
  seqnames = dmcs_delta_dt$chr,
  ranges = IRanges(start = dmcs_delta_dt$pos, end = dmcs_delta_dt$pos),
  pval = dmcs_delta_dt$pval
)

dmr_gr <- GRanges(
  seqnames = dmrs_char$chr,
  ranges = IRanges(start = dmrs_char$start, end = dmrs_char$end)
)

hits <- findOverlaps(dmc_gr, dmr_gr)

dmc_in_dmrs <- unique(queryHits(hits))

# Fraction of DMCs inside DMRs
frac_in_dmrs <- length(dmc_in_dmrs) / length(dmc_gr)

# Number of DMCs inside DMRs
length(dmc_in_dmrs)

# Isolated DMRs
dmr_with_dmc <- unique(subjectHits(hits))
isolated_dmrs <- setdiff(seq_along(dmrs$chr), dmr_with_dmc)
length(isolated_dmrs)

# Plot fraction
df <- data.frame(
  category = c("Inside DMR", "Outside DMR"),
  count = c(length(dmc_in_dmrs),
            length(dmc_gr) - length(dmc_in_dmrs))
)

ggplot(df, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Fraction of DMCs inside DMRs") +
  theme_bw()

ggsave("diagrams/fraction_dmcs_inside_dmrs.png",
       width = 7,
       height = 5,
       dpi = 300)

# DMR annotate promoters

if(!file.exists(paste0(outdatadir, "annotation_hg38_promoters.RDS"))) {
  annots_promoters <- c('hg38_genes_promoters')
  annotations_promoters <- build_annotations(genome = 'hg38', annotations = annots_promoters)
  saveRDS(annotations_promoters, paste0(outdatadir, "annotation_hg38_promoters.RDS")) } else {
  annotations_promoters <- readRDS(paste0(outdatadir, "annotation_hg38_promoters.RDS"))
  }

dmrs_gr <- GRanges(
  seqnames = dmrs_char$chr,
  ranges = IRanges(start = dmrs_char$start, end = dmrs_char$end),
  diff = dmrs_char$diff,
  nCG = dmrs_char$nCG
)

dmrs_annotated_promoters <- annotate_regions(
  regions = dmrs_gr,
  annotations = annotations_promoters,
  ignore.strand = TRUE,
  quiet = FALSE
)

dmrs_annot_promoters_dt <- as.data.table(dmrs_annotated_promoters)
dmrs_annot_promoters_dt[, annot_type := gsub("hg38_", "", annot.type)]
dmrs_annot_promoters_dt[, annot_type := gsub("genes_", "", annot_type)]
dmrs_annot_promoters_dt[, annot_type := gsub("cpg_", "", annot_type)]

dmrs_copy <- copy(dmrs_char)

# Take the first annotation for each DMC
dmrs_features_promoters <- dmrs_annot_promoters_dt[, .SD[1], by = .(seqnames, start)]
setnames(dmrs_features_promoters, c("seqnames"), c("chr"))

dmrs_final_promoters <- merge(dmrs_copy, 
                   dmrs_features_promoters[, .(chr, start, annot_type)],
                   by = c("chr", "start"), all.x = TRUE)

saveRDS(dmrs_final_promoters, paste0(outdatadir, "dmrs_annotated_promoters.RDS"))


# DMR Enrichment in promoters

setDT(dmrs_final_promoters)
feature_counts <- dmrs_final_promoters[, .N, by = annot_type]
feature_counts[, percentage := 100 * N / sum(N)]

dmrs_final_promoters[, feature_category := fcase(
  grepl("promoter", annot_type, ignore.case = TRUE), "Promoters",
  default = "Other"
)]

dmrs_promoters_feature_summary <- dmrs_final_promoters[, .(total_length = sum(length, na.rm = TRUE)), by = feature_category]
dmrs_promoters_feature_summary[, percentage := 100 * total_length / sum(total_length)]


background_promoters_dt <- as.data.table(annotations_promoters)
background_promoters_dt <- background_promoters_dt[seqnames == "chr20"]
background_promoters_dt[, annot_type := gsub("hg38_", "", type)]
background_promoters_dt[, annot_type := gsub("genes_", "", annot_type)]
background_promoters_dt[, annot_type := gsub("cpg_", "", annot_type)]


gr <- granges(BS_filtered)
cpgs_region_total_length <- max(end(gr)) - min(start(gr))
background_promoters_total_length <- sum(background_promoters_dt$width)
background_other_total_length <- cpgs_region_total_length - background_promoters_total_length

background_promoters_summary <- data.table(
  feature_category = c("Other", "Promoters"),
  total_length = c(
    background_other_total_length,
    background_promoters_total_length
  )
)

background_promoters_summary[, bg_percentage := 100 * total_length / sum(total_length)]

# Merge and calculate enrichment
promoters_enrichment_table <- merge(dmrs_promoters_feature_summary, 
                         background_promoters_summary[, .(feature_category, bg_percentage)],
                         by = "feature_category")

promoters_enrichment_table[, enrichment := percentage / bg_percentage]
promoters_enrichment_table[, log2_enrichment := log2(enrichment)]



# Diagrams

ggplot(dmrs_promoters_feature_summary, aes(x = reorder(feature_category, -percentage), 
                                   y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.5, size = 3.5) +
  labs(title = "DMR Total Length Distribution across Genomic Features",
       x = "Genomic Feature",
       y = "Percentage of DMCs (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = "diagrams/dmr_promoters_genomic_features.png",
  width = 8,
  height = 6,
  dpi = 300
)


ggplot(promoters_enrichment_table, aes(x = reorder(feature_category, enrichment), 
                                    y = enrichment)) +
  geom_bar(stat = "identity", aes(fill = enrichment > 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                    labels = c("Depleted", "Enriched")) +
  geom_text(aes(label = round(enrichment, 2)), hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "DMR Enrichment: Promoters vs Other",
       x = "Genomic Feature",
       y = "Enrichment (Observed / Expected)",
       fill = "") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  filename = "diagrams/dmr_promoters_enrichment.png",
  width = 8,
  height = 6,
  dpi = 300
)

# Features by direction

feature_by_direction <- dmrs_final_promoters[, .(total_length = sum(length, na.rm = TRUE)), by = .(feature_category, direction)]

feature_by_direction[, total := sum(total_length), by = direction]
feature_by_direction[, percentage := 100 * total_length / total]


ggplot(feature_by_direction, aes(x = feature_category, y = percentage, 
                                        fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Hypermethylated" = "#E64B35",
                               "Hypomethylated" = "#4DBBD5")) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = "Promoter Preferences: Hyper vs Hypo DMCs",
       x = "Genomic Feature",
       y = "Percentage of DMCs (%)",
       fill = "Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(
  filename = "diagrams/dmc_feature_by_direction.png",
  width = 8,
  height = 6,
  dpi = 300
)
