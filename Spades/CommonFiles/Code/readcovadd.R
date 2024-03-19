df<-read.csv("clean_contigs.stats.csv")
cov<-read.csv("read_coverage.tsv", sep = "\t")
df$ReadDepth<-as.numeric(cov[1,2])
write.csv(df, "clean_contigs.stats.csv", row.names = FALSE)