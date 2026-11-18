#Calculate Simpsons diversity in figure 2

library(vegan)
library(dplyr)

ranks <- list(
  Phylum = "combined_phyla_coverM.tsv",
  Class  = "combined_classes_coverM.tsv",
  Order  = "combined_orders_coverM.tsv"
)

calc_simpson <- function(file) {
  df <- read.table(file, sep = "\t", header = TRUE, row.names = 1)
  df[is.na(df)] <- 0
  D <- diversity(df, index = "simpson")
  data.frame(Sample = names(D), Simpson = D)
}

simpson_all <- bind_rows(lapply(ranks, calc_simpson), .id = "Rank")

write.csv(simpson_all, "Simpson_Diversity.csv", row.names = FALSE)
