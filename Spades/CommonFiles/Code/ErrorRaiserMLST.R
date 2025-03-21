local<-list.files(pattern = "localmlst.tsv")[1]
api<-list.files(pattern = "seqmlst.csv")[1]

local<-read.csv(local, sep = "\t", header = FALSE)
api<-read.csv(api)

if(is.na(api[1,1]) & ncol(local)==3){
  Sys.sleep(45)
  stop("Error: Something went wrong!")
}
  