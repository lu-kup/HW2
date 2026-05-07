library(data.table)
library(tidyverse)



sample_cols <- setdiff(colnames(coverage_dt), c("cpg_id", "chr", "position"))

coverage_long <- coverage_dt %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "coverage"
  )

coverage_long <- coverage_long %>%
  filter(!is.na(coverage))

head(coverage_long)

# Coverage distribution plot
p <- ggplot(coverage_long, aes(x = coverage)) +
  geom_histogram(
    binwidth = 2,
    fill = "#4C9BE8",
    color = "white"
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 300000)) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(
    title = "Coverage Distribution per Sample",
    x = "Coverage Depth",
    y = "Number of CpGs"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  "diagrams/coverage_distribution.png",
  plot = p,
  width = 12,
  height = 8
)


cpg_counts_thresh <- coverage_long %>%
  mutate(threshold = case_when(
    coverage >= 20 ~ ">=20x",
    coverage >= 10 ~ ">=10x",
    coverage >= 1 ~ ">=1x",
    TRUE ~ "All (>0)"
  )) %>%
  group_by(sample, threshold) %>%
  summarise(n_cpgs = n_distinct(cpg_id), .groups = "drop")

p2 <- ggplot(cpg_counts_thresh, aes(sample, n_cpgs, fill = threshold)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    breaks = function(x) seq(0, max(x, na.rm = TRUE), by = 50000)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Covered CpGs per Sample by Coverage Threshold",
    x = "Sample",
    y = "Number of CpGs",
    fill = "Threshold"
  )

ggsave(
  "diagrams/covered_cgs.png",
  plot = p2,
  width = 12,
  height = 8
)

summary1 <- summary(coverage_dt)

summary2 <- coverage_dt %>%
  summarise(
    number_covered_cpgs = sum(!is.na(value_column)),
    total_number_reads = sum(value_column, na.rm = TRUE)
  )

summary2 <- coverage_dt %>%
  summarise(
    count = n(),
    across(
      where(is.numeric),
      list(
        number_covered_cpgs = ~sum(!is.na(.x)),
        total_number_reads = ~sum(.x, na.rm = TRUE)
      )
    )
  )

summary2 <- summary2 %>% select(-(1:3))

write.csv(summary1, "outputs/summary1.csv", row.names = FALSE)
write.csv(summary2, "outputs/summary2.csv", row.names = FALSE)


# Correlation plot
mat <- methylation_dt %>%
  select(where(is.numeric)) %>%
  select(-position) %>%
  as.matrix()

cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

library(reshape2)
cor_long <- melt(cor_mat)

p3 <- ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Pearson Correlation Between Samples",
    x = "",
    y = "",
    fill = "Correlation"
  )

ggsave(
  "diagrams/correlation_plot.png",
  plot = p3,
  width = 12,
  height = 8
)

# PCA
mat_t <- t(mat) 
mat_clean <- mat_t[, apply(mat_t, 2, function(x) all(is.finite(x)))]
pca <- prcomp(mat_clean, center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df$Group <- ifelse(grepl("healthy", pca_df$Sample), "Healthy", "Tumor")
pca_df$Tissue <- ifelse(grepl("R2", pca_df$Sample), "R2", ifelse(grepl("R7", pca_df$Sample), "R7", NA))

p4 <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Tissue, label = Sample)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(R2 = 16, R7 = 17)) +
  theme_bw() +
  labs(
    title = "PCA of Methylation Profiles",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  "diagrams/pca.png",
  plot = p4,
  width = 12,
  height = 8
)
