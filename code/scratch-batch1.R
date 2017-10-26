library("dplyr")
x <- read.delim("../data/batch1_qc.txt", stringsAsFactors = FALSE)

out <- x %>%
  mutate(experiment = paste0("0", experiment),
         sample = paste(experiment, well, sep = "-"),
         tra1.60 = as.logical(tra1.60),
         individual.1 = paste0("NA", individual.1),
         individual.2 = paste0("NA", individual.2),
         individual.3 = paste0("NA", individual.3),
         individual.4 = paste0("NA", individual.4),
         ERCC = ifelse(is.na(ERCC), "Not added", ERCC),
         ERCC = factor(ERCC, levels = c("Not added", "1:100000", "1:50000"),
                             labels = c("Not added", "100x dilution", "50x dilution")),
         ERCC = as.character(ERCC)) %>%
  select(sample, experiment:index) %>%
  arrange(sample)

# Change individual NA18522 to NA18852
out <- out %>% mutate(individual.1 = ifelse(individual.1 == "NA18522",
                                            "NA18852", individual.1),
                      individual.2 = ifelse(individual.2 == "NA18522",
                                            "NA18852", individual.2),
                      individual.3 = ifelse(individual.3 == "NA18522",
                                            "NA18852", individual.3),
                      individual.4 = ifelse(individual.4 == "NA18522",
                                            "NA18852", individual.4))

# Change individual NA19201 to NA19092
out <- out %>% mutate(individual.1 = ifelse(individual.1 == "NA19201",
                                            "NA19092", individual.1),
                      individual.2 = ifelse(individual.2 == "NA19201",
                                            "NA19092", individual.2),
                      individual.3 = ifelse(individual.3 == "NA19201",
                                            "NA19092", individual.3),
                      individual.4 = ifelse(individual.4 == "NA19201",
                                            "NA19092", individual.4))

# Change individual NA18510 to NA18507
out <- out %>% mutate(individual.1 = ifelse(individual.1 == "NA18510",
                                            "NA18507", individual.1),
                      individual.2 = ifelse(individual.2 == "NA18510",
                                            "NA18507", individual.2),
                      individual.3 = ifelse(individual.3 == "NA18510",
                                            "NA18507", individual.3),
                      individual.4 = ifelse(individual.4 == "NA18510",
                                            "NA18507", individual.4))

str(out)

dir.create("../data/lab-info", showWarnings = FALSE)

for (e in unique(out$experiment)) {
  print(e)
  fname <- paste0("../data/lab-info/", e, ".txt")
  d <- out %>% filter(experiment == e)
  write.table(d, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
}
