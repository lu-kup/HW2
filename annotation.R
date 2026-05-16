library(pacman)
p_load(data.table, bsseq, ggplot2, DSS, data.table, annotatr, GenomicRanges, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, pheatmap)
options(scipen=999)
outdatadir <- "./outputs/"

# Annotation

# Build annotations - we set elements we need
if(!file.exists(paste0(outdatadir, "annotation_hg38.RDS"))) {
  annots <- c('hg38_genes_promoters')
  annotations <- build_annotations(genome = 'hg38', annotations = annots)
  saveRDS(annotations, paste0(outdatadir, "annotation_hg38.RDS")) } else {
  annotations <- readRDS(paste0(outdatadir, "annotation_hg38.RDS"))
  }
# Convert DMCs to GRanges as annotatr works with ranges.
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

# Convert and clean 
dmcs_annot_dt <- as.data.table(dmcs_annotated)
dmcs_annot_dt[, annot_type := gsub("hg38_", "", annot.type)]
dmcs_annot_dt[, annot_type := gsub("genes_", "", annot_type)]
dmcs_annot_dt[, annot_type := gsub("cpg_", "", annot_type)]

# Get gene information. Don't mix genome versions
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- genes(txdb)
# You can find any nearest thing using similar command
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

# Combine DMC data with nearest gene info
dmcs_with_genes <- copy(dmcs_delta_dt)
dmcs_with_genes[, dmc_index := 1:.N]

dmcs_with_genes <- merge(dmcs_with_genes, dmcs_nearest,
                         by = "dmc_index", all.x = TRUE)

# Add genomic feature annotations. Take the first annotation for each DMC
dmcs_features <- dmcs_annot_dt[, .SD[1], by = .(seqnames, start)]

setnames(dmcs_features, c("seqnames", "start"), c("chr", "pos"))

dmcs_final <- merge(dmcs_with_genes, 
                   dmcs_features[, .(chr, pos, annot_type)],
                   by = c("chr", "pos"), all.x = TRUE)

dmcs_final[, direction := ifelse(diff > 0, "Hypermethylated", "Hypomethylated")]

# Reorder columns
setcolorder(dmcs_final, c("chr", "pos", "diff", "pval", "fdr", "stat",
                         "direction", "annot_type", 
                         "gene_id", "nearest_gene_symbol", "nearest_gene_name",
                         "distance_to_gene"))

saveRDS(dmcs_annotated, paste0(outdatadir, "dmcs_annotated_granges.RDS"))
saveRDS(dmcs_final, paste0(outdatadir, "dmcs_final.RDS"))
