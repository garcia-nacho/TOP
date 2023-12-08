files.vec<-list.files()
tbp.files<-files.vec[grep("tbp_bulk.tar.gz", files.vec)]

#TBP
if(length(grep("dummy_tbp_bulk.tar.gz", files.vec))>0){
  file.remove(files.vec[grep("dummy_tbp_bulk.tar.gz", files.vec)])
}

files.vec<-list.files()
tbp.files<-files.vec[grep("tbp_bulk.tar.gz", files.vec)]

if(length(tbp.files)>0){
  tbpfileinfo<-as.data.frame(file.info(tbp.files))
  tbp.files<-tbp.files[order(tbpfileinfo$ctime, decreasing = TRUE)]
  
  dir.create("TB_Pipeline")
  dir.create("TB_Pipeline/COPY_TO_REPORTS")
  dir.create("TB_Pipeline/COPY_TO_TB_PIPELINE_DATABASE")
  
  for (i in 1:length(tbp.files)) {
    system(paste("tar -xvzf ",tbp.files[i]))
    tb_bulk_list<-list.files("tb_bulk",recursive = TRUE, full.names = TRUE)
    #The oldest one
    if(i == 1){
      
      file.rename(tb_bulk_list[grep("tb_bulk/TB_all", tb_bulk_list)], gsub("tb_bulk","TB_Pipeline",tb_bulk_list[grep("tb_bulk/TB_all", tb_bulk_list)]))
      file.rename("tb_bulk/Global_collection_tree.nwk","TB_Pipeline/Global_collection_tree.nwk" )
      file.rename("tb_bulk/snippy-core.log","TB_Pipeline/snippy-core.log" )
    }
    
    file.rename(tb_bulk_list[grep("tb_bulk/COPY_TO_REPORTS", tb_bulk_list)], gsub("tb_bulk","TB_Pipeline",tb_bulk_list[grep("tb_bulk/COPY_TO_REPORTS", tb_bulk_list)]))
    
    folder.to.create<-gsub("/.*","",gsub("tb_bulk/COPY_TO_TB_PIPELINE_DATABASE/","",tb_bulk_list[grep("tb_bulk/COPY_TO_TB_PIPELINE_DATABASE", tb_bulk_list)][1])) 
    dir.create(paste("TB_Pipeline/COPY_TO_TB_PIPELINE_DATABASE/", folder.to.create,sep = ""))
    file.rename(tb_bulk_list[grep("tb_bulk/COPY_TO_TB_PIPELINE_DATABASE", tb_bulk_list)], gsub("tb_bulk","TB_Pipeline",tb_bulk_list[grep("tb_bulk/COPY_TO_TB_PIPELINE_DATABASE", tb_bulk_list)]))
    
    R.utils::copyDirectory(paste("tb_bulk/",folder.to.create,sep = ""),paste("TB_Pipeline/",folder.to.create,sep = "") )
    R.utils::removeDirectory("tb_bulk", recursive=TRUE)
  }
  try(file.remove(tbp.files))
  
}

