library(ggpubr)
library(ggplot2)

depthfiles<-list.files(pattern = "depth.tsv")

for (i in 1:length(depthfiles)) {

  df<-read.csv(depthfiles[i],sep = "\t", header = FALSE)  
  
  plots<-list()
  contigs<-unique(df$V1)
  for (j in 1:length(contigs)) {
    plots[[j]]<-ggplot(df[which(df$V1 == contigs[j]),])+
      geom_line(aes(V2/1000,V3))+
      geom_smooth(aes(V2/1000,V3))+
      theme_minimal()+
      xlab("Position (Kbp)")+
      ylab("Depth")+
      ggtitle(contigs[j])+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  if(length(contigs)<=40){
    ggarrange(plotlist = plots, ncol = 4, nrow = 10)
    ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = 18)  
  }else if(length(contigs)>40 & length(contigs)<=44){
    ggarrange(plotlist = plots, ncol = 4, nrow = 11)
    ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = 20)
  }else if(length(contigs)>44 & length(contigs)<=48){
    ggarrange(plotlist = plots, ncol = 4, nrow = 12)
    ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = 22)
  }else if(length(contigs)>48 & length(contigs)<=52){
    ggarrange(plotlist = plots, ncol = 4, nrow = 13)
    ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = 24)
  }else if(length(contigs)>48 & length(contigs)<=52){
    ggarrange(plotlist = plots, ncol = 4, nrow = 14)
    ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = 26)
  }else{
      ggarrange(plotlist = plots, ncol = 4, nrow = round(length(contis)/4)+1)
      ggsave(gsub("_.*","_contigdepth.pdf", depthfiles[i]), width=12, height = (round(length(contis)/4)+1)*2)
    }
    
  }



  }


