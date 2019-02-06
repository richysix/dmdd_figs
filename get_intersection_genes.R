# get genes present in at least a supplied number of lines
Args <- commandArgs()
upset_wide <- read.delim(Args[6])
out_file <- Args[7]
threshold <- ifelse(is.na(Args[8]), length(upset_wide) - 1, as.integer(Args[8]) )
intersect_all <- upset_wide[ rowSums(upset_wide[,seq(2,length(upset_wide),1)]) >= threshold, ]
write.table(as.data.frame(intersect_all$gene), file = out_file, quote = FALSE,
row.names = FALSE, col.names = FALSE)
