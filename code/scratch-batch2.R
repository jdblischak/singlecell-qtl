library("dplyr")
x <- read.delim("../data/batch2_qc.txt", stringsAsFactors = FALSE)

out <- x %>%
  mutate(experiment = paste0("0", experiment),
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

dir.create("../data/lab-info", showWarnings = FALSE)

for (e in unique(out$experiment)) {
  print(e)
  fname <- paste0("../data/lab-info/", e, ".txt")
  d <- out %>% filter(experiment == e)
  write.table(d, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
}
