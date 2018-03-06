Args <- commandArgs()
upset_wide <- read.delim(Args[6])
intersect_all <- upset_wide[ rowSums(upset_wide[,seq(2,length(upset_wide),1)]) == length(upset_wide) - 1, ]
write.table(as.data.frame(intersect_all$gene), file = Args[7], quote = FALSE,
row.names = FALSE, col.names = FALSE)
