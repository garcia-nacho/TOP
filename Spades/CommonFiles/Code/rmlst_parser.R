library(jsonlite)


input<-list.files(pattern = "*_rmlst.json")
df<-fromJSON(input)

output<-df$taxon_prediction

#output$species<-df$fields$species
if(!is.null(df$fields$genus)){ 
output$genus<-df$fields$genus
}else{
output$genus<- gsub(" .*","",output$taxon)   
}
colnames(output)<-paste("rMLST_",colnames(output),sep="")

if(!is.null(df$fields$rST)){ 
output$rST<-df$fields$rST
}else{
output$rST<-NA
}

output$Sample<-gsub(".*/","",gsub("_.*","",input))

if(!is.null(df$exact_matches)){
 for (i in 1:length(df$exact_matches)) {
   dummy<-df$exact_matches[i][[1]] 
   AlleleID<-vector()
   if(!is.null(dummy$linked_data)){
   for (j in 1:nrow(dummy$linked_data)) {
     allele<-dummy$allele_id[j]
     sp<-dummy$linked_data$`rMLST genome database`$species[j][[1]]
     sp<-sp$value[which(sp$frequency==max(sp$frequency))[1]]
     AlleleID<-c(AlleleID,c(paste(sp, "/AlleleID:",allele,sep = "" )))
   }
   }else{
     AlleleID<-paste("ND/AlleleID:",dummy$allele_id,sep = "")
   }
  
   aid<-paste( nrow(dummy), "X ", names(df$exact_matches)[i],"/Sp:", paste(AlleleID, collapse = " | "),  sep = "")
   if(!exists("aid.out")){
     aid.out<-aid
   }else{
     aid.out<-c(aid.out,aid)
   }
 }
  if(exists("aid.out")){
    output$rMLST_scheme<-paste(aid.out, collapse = " | ")
  }
}

write.csv(output, gsub("_rmlst.json","_rmlst.csv", input),row.names = FALSE )