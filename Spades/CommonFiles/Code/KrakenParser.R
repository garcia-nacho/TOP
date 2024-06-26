library(data.table)
file<-list.files(pattern = ".*kraken.tsv")
df<-fread(file, sep = "\t")
df<-as.data.frame(df)
df<-aggregate(V2~V3, df, length)
colnames(df)<-c("Specie","Count")
df$Ratio<-df$Count/sum(df$Count)
df$Sample<-gsub("kraken.tsv", "", file)
df$Sample<-gsub("_.*","",df$Sample)
write.csv(df, gsub( "kraken.tsv","resultskraken.csv", file), row.names = FALSE)
file.remove(file)
