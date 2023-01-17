library(jsonlite)
library(ggplot2)

input<-list.files(pattern = "*_rmlst.json")



for (i in 1:length(input)) {
  
  df<-fromJSON(input[i])
  results<-as.data.frame(names(df$exact_matches))
  colnames(results)<-"Rmlst.alleles"
  results$Copies<-NA
  results$Sp<-NA
  results$Location<-NA
  results$AlleleID<-NA
  
  for (j in 1:nrow(results)) {
    results$Copies[j]<-length(  df$exact_matches[[j]]$allele_id) 
    results$AlleleID[j]<-paste(df$exact_matches[[j]]$allele_id, collapse = "/") 
    results$Location[j]<-paste(df$exact_matches[[j]]$contig, collapse = "/")
    sp<-vector()
    for (l in 1:length(  df$exact_matches[[j]]$allele_id) ) {
      if(!is.null(df$exact_matches[[j]]$linked_data)){
      sp.df<-df$exact_matches[[j]]$linked_data$`rMLST genome database`$species[[l]]
      sp<-c(sp,sp.df$value[which(sp.df$frequency==max(sp.df$frequency))])}
    }
    if(length(sp)>0){
    results$Sp[j]<-paste(sp, collapse = "/")}
  }
  results$Sample<-gsub("_.*","",input[i])

  if(!exists("results.out")){
    results.out<-results
  }else{
    results.out<-rbind(results.out,results)
  }
    
}




ggplot(results.out)+
  geom_point(aes(Rmlst.alleles, Copies,col=Sp),alpha=1)+
  theme_minimal()+
  scale_color_manual(values = rainbow(length(unique(results.out$Sp))))+
  ylab("Copy number")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Sample)


results.out$Sp[grep("/",results.out$Sp)]<-"Mix Sp"

ggplot(results.out)+
  geom_point(aes(Rmlst.alleles, Copies,col=Sp),alpha=1)+
  theme_minimal()+
  scale_color_manual(values = rainbow(length(unique(results.out$Sp))))+
  ylab("Copy number")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Sample)
