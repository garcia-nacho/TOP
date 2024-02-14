
compressed<-list.files(pattern = ".tar.gz")

if(length(grep("dummy_tbp", compressed))>0){
 file.remove(compressed[grep("dummy_tbp", compressed)])
 compressed<-list.files(pattern = ".tar.gz")
}

if(length(compressed)>0){
  for (i in 1:length(compressed)) {
    system(paste("tar -xf",compressed[i]))
  }
  file.remove(compressed)
}


samples<-list.dirs(recursive = FALSE)

if(length(samples)>0){
samples.df <- as.data.frame(samples)
samples.df$dirs<-gsub("./","",samples.df$samples)
samples.df$samples<-gsub("_.*","", samples.df$dirs)
samples.df$mash<-NA
if(exists("mykrobe")) rm(mykrobe)
for (i in 1:nrow(samples.df)) {
  my.dum<-read.csv(paste(samples.df$dirs[i],"/mykrobe_output.csv",sep = ""))
  my.dum<-my.dum[1,]
  if(!exists("mykrobe")){
    mykrobe<-my.dum
  }else{
    mykrobe<-rbind(mykrobe,my.dum)
  }
  
}

mykrobe<-mykrobe[,c("sample", "phylo_group", "species", "lineage" )]

for (i in 1:nrow(samples.df)) {
  error<-list.files(path = samples.df$dirs[i], pattern = "Mashclassificationproblem|Mashothermycobacterium")
  if(length(error)==1){
    error.df<-read.csv(paste(samples.df$dirs[i],"/",error,sep = ""), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    samples.df$mash[i]<-as.character(error.df[1,1])
  }else{
    samples.df$mash[i]<-"MTBC"
  }

}


mykrobe$sample<-gsub("_.*","", mykrobe$sample)

samples<-merge(samples.df, mykrobe, by.x = "samples", by.y = "sample")
samples$Excluded<-"YES"
if(length(which(samples$phylo_group=="Mycobacterium_tuberculosis_complex"))>0){
  samples$Excluded[which(samples$phylo_group=="Mycobacterium_tuberculosis_complex")]<-"NO"
}

to.remove<-samples$dirs[which(samples$Excluded=="YES")]
if(length(to.remove)>0){
  for (i in 1:length(to.remove)) {
    system(paste("rm -rf ",to.remove[i]))
  }
}
samples$dirs<-NULL
write.table(samples, "mashmykrobe.csv", row.names = FALSE, sep = "\t")

}else{
  write.table("NonTB_in_the_run","mashmykrobe.csv", row.names = FALSE, sep = "\t")
}
