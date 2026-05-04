
install.packages("ggplot2")
library(ggplot2)


cpg <- read.table("outputs/SRR11647649_chr20_CpG.bedGraph", header = FALSE,
  comment.char = "t")

# Assign column names (standard MethylDackel format)
colnames(cpg) <- c("chr", "start", "end", "meth_percent", "meth", "unmeth")

# -------------------------------
# 2. Compute coverage per CpG
# -------------------------------
cpg$coverage <- cpg$meth + cpg$unmeth

# -------------------------------
# 3. Total reads covering CpGs
# -------------------------------
total_reads_cpg <- sum(cpg$coverage)

cat("Total reads covering CpGs:", total_reads_cpg, "\n")

# -------------------------------
# 4. Coverage distribution
# -------------------------------

# Summary statistics
cat("\nCoverage summary:\n")
print(summary(cpg$coverage))

# Histogram (base R)
hist(cpg$coverage,
     breaks = 50,
     main = "Coverage distribution per CpG",
     xlab = "Coverage")

# Histogram (ggplot2)
ggplot(cpg, aes(x = coverage)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = "CpG Coverage Distribution",
       x = "Coverage",
       y = "Count")


write.table(cpg,
            file = "CpG_with_coverage.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save histogram
ggsave("CpG_coverage_distribution.png")

# ===============================
# End of script
# ===============================