require(dplyr)
require(data.table)

# read file path
all_paths <-
  list.files(path = "survival_results/",
             pattern = "*.csv",
             full.names = TRUE)

# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = ",",
         encoding = "UTF-8")

# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
# change column name
names(all_result)[3] <- "File.Path"
write.csv(all_result, file = "Combined.csv")

# Create relative abundance matrix
df1 <- read.csv("Combined.csv")
df2 <- with(df1, tapply(relative_abundance, list(clade_name, SampleID), sum))
df2[is.na(df2)] <- 0
write.csv(df2, file = "microbiomeWGS_matrix_a_matrix.csv")