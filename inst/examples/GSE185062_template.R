# GSE185062-style analysis template using detectPanel
# Adjust file paths and metadata column names to your local files.

library(detectPanel)

counts_df <- read.csv("GSE185062_miRNA_counts_matrix.csv", check.names = FALSE)
rownames(counts_df) <- counts_df$Gene
counts_df$Gene <- NULL
counts <- as.matrix(counts_df)

meta <- read.csv("GSE185062_sample_info.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$GSM

keep <- intersect(colnames(counts), rownames(meta))
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
meta <- meta[meta$group %in% c("Control", "NAFL"), , drop = FALSE]
counts <- counts[, rownames(meta), drop = FALSE]

result <- discover_panel(
  counts = counts,
  metadata = meta,
  outcome = "group",
  positive = "NAFL",
  exclude = "miR-3687",
  exclude_regex = "^(CTRL_|HK_)",
  panel_size = 3,
  candidate_n = 12,
  detection_threshold = 10,
  min_mean = 100,
  min_median = 20,
  min_detection = 0.80,
  min_group_detection = 0.50,
  min_auc = 0.80,
  outer_v = 5,
  outer_repeats = 5,
  inner_v = 5,
  inner_repeats = 10,
  allow_fallback = FALSE,
  seed = 123
)

print(result)
print(summary(result))
plot(result, type = "roc")
plot(result, type = "feature_frequency")
export_results(result, "GSE185062_detectPanel_results")
