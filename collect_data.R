sample_info <- data.frame(
  sample_id = c("T1_healthy_R2", "T1_healthy_R7",
                "T1_tumor_R2", "T1_tumor_R7"),
  tissue = rep(c("T1", "T1", "T1", "T1"), 1),
  condition = rep(c("healthy", "healthy", 
                    "tumor", "tumor"), 1),
  file_path = c("SRR11647649_chr20_CpG.bedGraph", "SRR11647651_chr20_CpG.bedGraph", 
                "SRR11647659_chr20_CpG.bedGraph", "SRR11647661_chr20_CpG.bedGraph"),
  stringsAsFactors = FALSE
)

sample_info$color <- ifelse(sample_info$tissue == "T1", 
                             ifelse(sample_info$condition == "healthy", "#4DBBD5", "#3C8DAD"),
                             ifelse(sample_info$condition == "healthy", "#E64B35", "#C43828"))
sample_info$group <- gsub("_R1|_R2|_R3|_R7", "", sample_info$sample_id)

cat("Sample Information:\n")
print(sample_info[, c("sample_id", "tissue", "condition")])
cat("\n")

library(data.table)
outdatadir <- "outputs/"

if(!((file.exists(paste0(outdatadir, "coverage_dt.RDS"))) & (file.exists(paste0(outdatadir, "methylation_dt.RDS"))) & (file.exists(paste0(outdatadir, "sample_info.RDS"))))) {
  idir <- "outputs/"
  data_list <- list()
  for(i in 1:nrow(sample_info)) {  
    tmp <- read.table(paste0(idir, sample_info$file_path[i]), 
                      header = FALSE,
                      skip = 1,
                      stringsAsFactors = FALSE)                 

    setnames(tmp, c("chr", "position", "end_position",
                    "rounded_methylation",
                    "methylated_read_count", 
                    "unmethylated_read_count"))

    tmp <- as.data.table(tmp)

    data_list[[sample_info$sample_id[i]]] <- tmp
  }

  # We will create two matrices: one for methylation levels and one for coverage


  for(i in 1:length(data_list)) {
    data_list[[i]][, cpg_id := paste(chr, position, sep = ":")]
  }
  all_cpgs <- data_list[[1]]$cpg_id
  for(i in 2:length(data_list)) {
    all_cpgs <- union(all_cpgs, data_list[[i]]$cpg_id)
  }

  all_cpgs_df <- data.table(cpg_id = all_cpgs)
  all_cpgs_df[, c("chr", "position") := tstrsplit(cpg_id, ":", fixed = TRUE)]
  all_cpgs_df[, position := as.integer(position)]

  reference_dt <- all_cpgs_df
  setkey(reference_dt, cpg_id)

  methylation_dt <- copy(reference_dt)
  for(i in 1:length(data_list)) {
    sample_name <- names(data_list)[i]
    
    sample_dt <- data_list[[i]][, .(cpg_id, methylated_read_count, unmethylated_read_count)]
    setkey(sample_dt, cpg_id)
    
    sample_dt[, beta := ifelse(methylated_read_count + unmethylated_read_count > 0, 
                                methylated_read_count / (methylated_read_count + unmethylated_read_count), 
                                NA_real_)]
    methylation_dt[sample_dt, (sample_name) := i.beta, on = "cpg_id"]
  }

  coverage_dt <- copy(reference_dt)
  for(i in 1:length(data_list)) {
    sample_name <- names(data_list)[i]
    
    sample_dt <- data_list[[i]][, .(cpg_id, methylated_read_count, unmethylated_read_count)]
    setkey(sample_dt, cpg_id)

    sample_dt[, total_read_count := methylated_read_count + unmethylated_read_count]
    coverage_dt[sample_dt, (sample_name) := i.total_read_count, on = "cpg_id"]
  }

  saveRDS(coverage_dt, paste0(outdatadir, "coverage_dt.RDS"))
  saveRDS(methylation_dt, paste0(outdatadir, "methylation_dt.RDS"))
  saveRDS(sample_info, paste0(outdatadir, "sample_info.RDS")) } else {

  coverage_dt <- readRDS(paste0(outdatadir, "coverage_dt.RDS"))
  methylation_dt <- readRDS(paste0(outdatadir, "methylation_dt.RDS"))
  sample_info <- readRDS(paste0(outdatadir, "sample_info.RDS"))    
  }
