library("reshape2")
Args            <- commandArgs()
data <- read.delim(Args[6])
count_col_names <- colnames(data)[ grepl("normalised\\.count$", colnames(data)) ]
# subset data to counts
count_data <- as.matrix(data[, count_col_names])
# centre and scale by row
count_data_scaled <- t(scale(t(count_data)))
# add back gene name and gene id columns
scaled_data <- cbind(data[, c("Gene.ID", "Name")], count_data_scaled)
colnames(scaled_data) <- gsub("\\.normalised\\.count$", "", colnames(scaled_data))
count_col_names <- gsub("\\.normalised\\.count$", "", count_col_names)
# melt
data_m <- melt(scaled_data, id.vars = c("Gene.ID", "Name"),
                measure.vars = count_col_names,
                variable.name = "Sample",
                value.name = "Normalised.Counts.Scaled")
write.table(data_m, file = Args[7], quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

