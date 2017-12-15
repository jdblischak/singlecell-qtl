library("dplyr")
x <- read.delim("../data/batch3_qc.txt", stringsAsFactors = FALSE)

out <- x %>%
  mutate(experiment = sprintf("%08d", experiment),
         sample = paste(experiment, well, sep = "-"),
         tra1.60 = as.logical(tra1.60),
         fly = 0,
         worm = 0,
         ERCC = ifelse(is.na(ERCC), "Not added", ERCC),
         ERCC = factor(ERCC, levels = c("Not added", "1:100000", "1:50000"),
                             labels = c("Not added", "100x dilution", "50x dilution")),
         ERCC = as.character(ERCC)) %>%
  select(sample, experiment:individual.4, fly, worm, ERCC, index) %>%
  arrange(sample)

# Confirm the leading zero was added properly now that we have dates from Oct
stopifnot(out$experiment %>% sort %>% unique %>% nchar == 8)

dir.create("../data/lab-info", showWarnings = FALSE)

for (e in unique(out$experiment)) {
  print(e)
  fname <- paste0("../data/lab-info/", e, ".txt")
  d <- out %>% filter(experiment == e)
  write.table(d, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
}

# For pasting into config.yaml:
# out$experiment %>% sort %>% unique %>% paste(collapse = '", "') %>% cat
