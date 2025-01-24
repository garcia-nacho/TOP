setwd("/media/nacho/Data/temp/toptest/BPE/TOPresults/stitching/PanGenomeBlastDB")

fil<-list.files(pattern = "Stitched", full.names = TRUE)
plotlist<-list()
for (i in 1:length(fil)) {
  dum<-read.csv(fil[i],sep = "\t")
  
  ice<-as.data.frame(table(dum$SubjectID))
  to.ignore<-as.character(ice$Var1[which(ice$Freq>150)])
  plotlist<-c(plotlist, list(ggplot(dum)+
                               geom_segment(aes(x=StartRef,xend=EndRef,y="Ref",yend="Ref"))+
                               geom_segment(aes(x=CorrectedStart,xend=CorrectedEnd,y="Query",yend="Query"))+
                               geom_segment( data=dum[-which(dum$SubjectID %in% to.ignore),],aes(x=CorrectedStart, xend= StartRef, y="Query", yend="Ref"), col="red", alpha=0.1)+
                               geom_segment( data=dum[-which(dum$SubjectID %in% to.ignore),], aes(x=CorrectedEnd, xend= EndRef, y="Query", yend="Ref"), col="red", alpha=0.1)+
                               theme_minimal()+
                               xlab("Reference Position")+
                               ylab("Genome") +
                               ggtitle(gsub("_clean.*","",gsub(".*/", "", fil[i])))))
  
  

}


ggarrange(plotlist = plotlist)



starts<-dum$StartRef[which(dum$SubjectID %in% to.ignore)]
starts<-starts[-which(is.na(starts))]
starts<-starts[-which(duplicated(starts))]

plot(density(starts))


ggplotly(ggplot(dum)+
           geom_segment(aes(x=StartRef,xend=EndRef,y="Ref",yend="Ref"))+
           geom_segment(aes(x=CorrectedStart,xend=CorrectedEnd,y="Query",yend="Query"))+
           geom_segment( aes(x=CorrectedStart, xend= StartRef, y="Query", yend="Ref"), col="red", alpha=0.1)+
           geom_segment( aes(x=CorrectedEnd, xend= EndRef, y="Query", yend="Ref"), col="red", alpha=0.1)+
           theme_minimal()+
           xlab("Reference Position")+
           ylab("Genome"))



headers<-list.files("/media/nacho/Data/OnGoingProjects/Pertusis/fastq", pattern = "headersR", full.names = TRUE, recursive = TRUE)

samples<-unique(gsub("_headersR.*", "", headers))

if(exists("out"))rm(out)

for (i in 1:length(samples)){

  r1.full<-read.csv(paste(samples[i],"_headersR1.tsv",sep = ""),sep = "\t", header = FALSE)
  r2.full<-read.csv(paste(samples[i],"_headersR2.tsv",sep = ""),sep = "\t", header = FALSE)
  
  genes<-unique(r1.full$V3)
  
  for (j in 1:length(genes)){
    r1<-r1.full[which(r1.full$V3==genes[j]),]
    r2<-r2.full[which(r2.full$V3==genes[j]),]
    r1frac<-length(which(r1$V1 %in% r2$V1))/nrow(r1)
    r2frac<- length(which(r2$V1 %in% r1$V1))/nrow(r2)
    temp<-as.data.frame(t(c(gsub(".*/","",samples[i]), 1-r1frac, 1-r2frac, nrow(r1), nrow(r2))))
    colnames(temp)<-c("Sample", "R1Unmapped","R2Unmapped","R1Count","R2Count")
    temp$Gene<-genes[j]
    if(!exists("out")){
      out<-temp
    }else{
      out<-rbind(out,temp)
    }
  }
}

out$Sample<-gsub("_.*","",out$Sample)

out$R1Unmapped<-as.numeric(out$R1Unmapped)

ggplot(out[which(out$Gene=="prn"),])+
  geom_bar(aes(Sample,R1Unmapped),stat = "identity")+
  theme_minimal()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df<-read.csv("/media/nacho/Data/OnGoingProjects/Pertusis/BlastRoI.csv")



# IS418 ---------------------------------------------------------------------

headers<-list.files("/media/nacho/Data/OnGoingProjects/Pertusis/fastq", pattern = "headersR*", full.names = TRUE, recursive = TRUE)

samples<-unique(gsub("_.*", "", headers))

if(exists("out"))rm(out)

pb<-txtProgressBar(max = length(samples))
for (i in 1:length(samples)){
  setTxtProgressBar(pb,i)
  dummyfiles<-list.files(paste(gsub( "fastq/.*","fastq/", samples[i]), gsub( ".*/","", samples[i]) ,sep = ""),pattern = "tsv",full.names = TRUE)
  
  R1s<-dummyfiles[grep("headersR1",dummyfiles)]
  ctrl<-dummyfiles[grep("IS418_headersR2.tsv",dummyfiles)]
  R1s<-R1s[-grep("IS418", R1s)] 
  R1ctrl<-read.csv(ctrl,sep = "\t", header = FALSE)
  
  for (j in 1:length(R1s)) {
    subject<-read.csv(R1s[j], sep = "\t", header = FALSE)
    frac<-length(which(subject$V1 %in% R1ctrl$V1))/nrow(subject)
    temp<-as.data.frame(t(c(frac, gsub(".*/","",R1s[j])) ))
    colnames(temp)<-c("Fraction","Sample")
    temp$Read<-"R1"
    if(!exists("out")){
      out<-temp
    }else{
      out<-rbind(out,temp)
    }
  
  }
  
  R1s<-dummyfiles[grep("headersR2",dummyfiles)]
  ctrl<-dummyfiles[grep("IS418_headersR1.tsv",dummyfiles)]
  R1s<-R1s[-grep("IS418", R1s)] 
  R1ctrl<-read.csv(ctrl,sep = "\t", header = FALSE)
  
  for (j in 1:length(R1s)) {
    subject<-read.csv(R1s[j], sep = "\t", header = FALSE)
    frac<-length(which(subject$V1 %in% R1ctrl$V1))/nrow(subject)
    temp<-as.data.frame(t(c(frac, gsub(".*/","",R1s[j])) ))
    colnames(temp)<-c("Fraction","Sample")
    temp$Read<-"R2"
    if(!exists("out")){
      out<-temp
    }else{
      out<-rbind(out,temp)
    }
    
  }
  
}


out$Sample<-gsub("_headers.*","",out$Sample)
out$Gene<-gsub("_Area",".Area",out$Sample)
out$Gene<- gsub(".*_","",out$Gene) 
out$Sample<-gsub("_.*","",out$Sample)

out$Fraction<-as.numeric(out$Fraction)
outagg<-aggregate(Fraction ~ Sample + Gene, out, mean)

ggplot(out)+
  geom_bar(aes(Sample,Fraction),stat = "identity")+
  theme_minimal()+
  facet_wrap(~Gene+Read)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(outagg)+
  geom_bar(aes(Sample,Fraction),stat = "identity")+
  theme_minimal()+
  facet_wrap(~Gene)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

nzgb<-genebankreader("/media/nacho/Data/DockerImages/TOP_dev/BPProfiler/Refs/NZ_CP026996.gb")

#1083517-1084467-1085772-1087193
#
sq<-read.fasta("/media/nacho/Data/OnGoingProjects/Pertusis/fastq/NZ_CP026996.fasta")
sq<-sq[[1]][1083517:1087193]
write.fasta(sq, "/media/nacho/Data/OnGoingProjects/Pertusis/PertactinInsertion.fasta", names="prn_IS418")

inserted<-sq[[1]][1084467:1085772]
write.fasta(inserted, "/media/nacho/Data/OnGoingProjects/Pertusis/Insertion.fasta", names="Insertion")

dp<-list.files("/media/nacho/Data/OnGoingProjects/Pertusis/fastq", pattern = "_NZ.depth.tsv", full.names = TRUE, recursive = TRUE)
outlist<-list()
for (i in 1:length(dp)) {
  dum<-read.csv(dp[i],sep = "\t", header = FALSE)
  outlist<-c(outlist,list(ggplot(dum[which(dum$V2>1083507 & dum$V2 < 1083627),])+
                            geom_line(aes(V2,V3))+
                            #geom_vline(xintercept = c(1083517,1084467,1085772,1087193),col="red")+
                            theme_minimal()+
                            ylim(0,max(dum$V3[which(dum$V2>1082000 & dum$V2 < 1088193)])+10)+
                            ggtitle(gsub("_NZ.de.*","",gsub(".*/","",dp[i])))))
  
    
}

ggarrange(plotlist = outlist)



a<-read.fasta("/media/nacho/Data/DockerImages/TOP_dev/BPProfiler/Refs/TohamaI.fasta")
a<-a[[1]][c(1097091:1101823)]
write.fasta(a, "/media/nacho/Data/DockerImages/TOP_dev/BPProfiler/Refs/IndividualSeqs/prn_LargeArea.fasta",names="prn_LargeArea")
