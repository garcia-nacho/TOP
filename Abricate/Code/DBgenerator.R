library(seqinr)

fastas<-read.fasta("/media/nacho/Data/DockerImages/TOP_dev/Abricate/vfdb/sequences")

names<-system("grep '^>' /media/nacho/Data/DockerImages/TOP_dev/Abricate/vfdb/sequences", intern = TRUE)

#Format
#">vfdb2~~~VFG037170(gb|WP_001081754)~~~(plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii 1656-2]"


names<-gsub(">", "", names)

counts<-vector()
for (i in 1:length(names)) {
  counts<-c(counts,length(gregexpr("\\) \\(", names[i])[[1]]) )
}

ids<-gsub("\\) \\(.*", ")", names )
ids<-gsub("\\) .*", ")",ids)
genes<-sub(".*\\) \\(", "(", names )

index<-which(counts==2)

for (i in 1:length(index)) {
  split_position <- regexpr("\\) ", names[index[i]])[1]
  
  # Split the string at that position
  ids[index[i]] <- substr(names[index[i]], 1, split_position)
  genes[index[i]] <- substr(names[index[i]], split_position + 2, nchar(names[index[i]]))  # +2 to remove the ") "
  
}


index<-which(counts==3)

for (i in 1:length(index)) {
  split_position <- regexpr("\\) ", names[index[i]])[1]
  
  # Split the string at that position
  ids[index[i]] <- substr(names[index[i]], 1, split_position)
  genes[index[i]] <- substr(names[index[i]], split_position + 2, nchar(names[index[i]]))  # +2 to remove the ") "
  
}

index<-which(counts==4)
for (i in 1:length(index)) {
  split_position <- regexpr("\\) ", names[index[i]])[1]
  
  # Split the string at that position
  ids[index[i]] <- substr(names[index[i]], 1, split_position)
  genes[index[i]] <- substr(names[index[i]], split_position + 2, nchar(names[index[i]]))  # +2 to remove the ") "
  
}



index<-which(counts==5)
for (i in 1:length(index)) {
  split_position <- regexpr("\\) ", names[index[i]])[1]
  
  # Split the string at that position
  ids[index[i]] <- substr(names[index[i]], 1, split_position)
  genes[index[i]] <- substr(names[index[i]], split_position + 2, nchar(names[index[i]]))  # +2 to remove the ") "
  
}


index<-c(1:length(genes))[-grep("^\\(", genes)]

for (i in 1:length(index)) {
  split_position <- regexpr("\\) ", names[index[i]])[1]
  
  # Split the string at that position
  ids[index[i]] <- substr(names[index[i]], 1, split_position)
  temp<-paste("(",gsub(" .*", "",ids[index[i]]),")", sep = "")
  genes[index[i]] <- paste(temp, substr(names[index[i]], split_position + 2, nchar(names[index[i]])))  # +2 to remove the ") "
  
}

seqnames<-paste("vfdb2", ids, genes, sep = "~~~")
write.fasta(fastas, "/media/nacho/Data/DockerImages/TOP_dev/Abricate/vfdb/sequences", names=seqnames)
