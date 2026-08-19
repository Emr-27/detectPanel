make_synthetic_panel_data <- function(n = 30L, p = 10L, seed = 1L) {
  set.seed(seed)
  y <- rep(c("Control", "Case"), length.out = n)
  counts <- matrix(
    rpois(p * n, lambda = 60),
    nrow = p,
    dimnames = list(paste0("m", seq_len(p)), paste0("s", seq_len(n)))
  )
  counts[1:3, y == "Case"] <- counts[1:3, y == "Case"] + 70
  metadata <- data.frame(group = y, row.names = colnames(counts))
  list(counts = counts, metadata = metadata, y = y)
}
