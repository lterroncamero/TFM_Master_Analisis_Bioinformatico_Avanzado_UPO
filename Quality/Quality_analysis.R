#!/usr/bin/env Rscript

format_enrichment<-function(ego,de_data){
  data<-data.frame(ego)
  colnames(data)<-gsub("Count","Total",colnames(data))
  genes_down<-de_data[de_data$FDR<=0.05 & de_data$logFC<0,"Gene"]
  genes_up<-de_data[de_data$FDR<=0.05 & de_data$logFC>0,"Gene"]
  data$GenesID_up<-NA
  data$GenesID_dn<-NA
  data$Total_up<-NA
  data$Total_dn<-NA
  
  for(n in 1:nrow(data)){
    up<-genes_up[genes_up %in% unlist(strsplit(data[n,"geneID"],"/"))]
    dn<-genes_down[genes_down %in% unlist(strsplit(data[n,"geneID"],"/"))]
    data[n,"GenesID_up"]<-paste(up,collapse = "/")
    data[n,"GenesID_dn"]<-paste(dn,collapse = "/")
    data[n,"Total_up"]<-length(up)
    data[n,"Total_dn"]<-length(dn)
  }
  return(data[,c(1:8,10,11,9,12,13)])
}
my_enrichment<-function(de_data,FDR=0.05,FA_label,cutoff=0.05,organism="human",FC=0){
  de_data<-probeset.list
  FDR<-0.05
  FA_label<-"LH7vsLH2"
  cutoff<-0.05
  organism<-"human"
  FC<-0
  
  dbName<-NULL
  keggDB<-NULL
  if(organism == "human"){
    library("org.Hs.eg.db")
    dbName<-org.Hs.eg.db
    keggDB<-"hsa"
  }else if (organism == "mouse"){
    library("org.Mm.eg.db")
    dbName<-org.Mm.eg.db
    keggDB<-"mmu"
    
  }else{
    stop("At the moment only human and mouse are supported")
  }
  
  submod = read.delim("/mnt/beegfs/Bioinformatica/EXT/08-EXT-20-Saltmae-endometrio/Analisis/arrays/WGCNA/WGCNA_modules_array_Quality.tsv")
  submod$gene_id
  turquoise<-submod[submod$colors == 'turquoise',]
  gene<-turquoise$gene_id
  
  #gene<-de_data[de_data$adj.P.Val<=FDR & abs(de_data$logFC)>=FC,"Group.1"]
  
  eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=dbName)
  all_eg = bitr(submod$gene_id, fromType="SYMBOL", toType="ENTREZID", OrgDb=dbName)
  
  new_path=paste(path,"Enrichment/",sep="")
  new_path<-"/mnt/beegfs/Bioinformatica/EXT/08-EXT-20-Saltmae-endometrio/Analisis/arrays/WGCNA/Enrichment/"
  system(paste("mkdir -p ", new_path,sep=""))
  ego_BP <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  ego_MF <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  
  ego_CC <- enrichGO(gene          = as.vector(eg$ENTREZID),
                     universe      = as.vector(all_eg$ENTREZID),
                     OrgDb         = dbName,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = cutoff,
                     qvalueCutoff  = cutoff,
                     readable      = TRUE)
  FA_label<-"turquoise"
  if(nrow(as.data.frame(ego_BP))>0){
    write.table(format_enrichment(ego_BP,de_data),file = paste(new_path,"Enrichment_GO_BP_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_BP_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_BP)
    y$labels$title<-"GO terms from BP\n"
    plot(y)
    
    x<-dotplot(ego_BP)
    x$labels$title<-"GO terms from BP\n"
    plot(x)
    
    #nose<-ego_BP
    #class(nose)<-NULL
    #z<-cnetplot(ego_BP,categorySize = "pvalue",main="GO terms from BP\n",vertex.label.cex=0.5)
    #plot(z)
    
    plotGOgraph(ego_BP, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste(new_path,"GO_BP_dotplot",FA_label,".png",sep=""))
    plot(x)
    dev.off()
  }
  if(nrow(as.data.frame(ego_MF))>0){
    
    write.table(format_enrichment(ego_MF,de_data),file = paste(new_path,"Enrichment_GO_MF_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_MF_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_MF)
    y$labels$title<-"GO terms from MF\n"
    plot(y)
    
    x<-dotplot(ego_MF)
    x$labels$title<-"GO terms from MF\n"
    plot(x)
    
    #z<-cnetplot(ego_MF,categorySize = "pvalue",main="GO terms from MF\n",vertex.label.cex=0.5)
    #plot(z)
    plotGOgraph(ego_MF, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste(new_path,"GO_MF_dotplot",FA_label,".png",sep=""))
    plot(x)
    dev.off()
  }
  if(nrow(as.data.frame(ego_CC))>0){
    write.table(format_enrichment(ego_CC,de_data),file = paste(new_path,"Enrichment_GO_CC_",FA_label,".xls",sep=""),sep="\t",row.names = F)
    pdf(file = paste(new_path,"Enrichment_GO_CC_",FA_label,".pdf",sep=""),paper="a4")
    y<-barplot(ego_CC)
    y$labels$title<-"GO terms from CC\n"
    plot(y)
    
    x<-dotplot(ego_CC)
    x$labels$title<-"GO terms from CC\n"
    plot(x)
    
    #z<-cnetplot(ego_CC,categorySize = "pvalue",main="GO terms from CC\n",vertex.label.cex=0.3,node.color="red",edge.width=0.2)
    #plot(z)
    plotGOgraph(ego_CC, firstSigNodes = 5,  sigForAll = TRUE, useFullNames = TRUE)
    dev.off()
    
    png(file = paste(new_path,"GO_CC_dotplot",FA_label,".png",sep=""))
    plot(x)
    dev.off()
  }
  
  KEGG_enrich <- enrichKEGG(
    gene= as.vector(eg$ENTREZID),
    keyType = "ncbi-geneid",
    organism     = keggDB,
    pvalueCutoff = cutoff,
    use_internal_data = F
  )
  if(is.null(KEGG_enrich) == F){
    if(nrow(KEGG_enrich)>0){
      write.table(KEGG_enrich,file = paste(new_path,"Enrichment_KEGG_",FA_label,".xls",sep=""),sep="\t",row.names = F)
      
      result<-NA
      for(cont in 1:length(KEGG_enrich$ID)){
        cat(paste("For pathway named: ", KEGG_enrich$ID[cont]))
        result[cont]<-browseKEGG(KEGG_enrich, KEGG_enrich$ID[cont])
        #print(result[cont])
        Sys.sleep(2)
        file<-readLines(result[cont])
        #cat(file)
        pattern='<.+ src=\"(.+)\" srcset='
        #cat(pattern)
        path<-gsub(pattern,"\\1", grep("pathwayimage",file,value=TRUE))
        #cat(path)
        result[cont]<-strsplit(path," ")[[1]][1]
        #cat(result[cont])
        result[cont]<-gsub("\\\"","",result[cont])
        result[cont]<-gsub("\\t","",result[cont])
        result[cont]<-gsub("\\.png.*",".png",result[cont])
        #cat(result[cont])
        my_url<-paste("http://www.kegg.jp",result[cont],sep="")
        #cat(my_url)
        print(paste("Downloading ", my_url))
        curl_download(my_url,paste(KEGG_enrich$ID[cont],".png",sep=""),mode="wb")
      }
      
      files <- list.files(path=".", pattern="*.png", all.files=T, full.names=T)
      
      pdf(file = paste(new_path,"Enrichment_KEGG_",FA_label,".pdf",sep=""),paper = "a4")
      for(i in 1:length(files)) {
        img <- readPNG(files[i])
        plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n",axes=F,xlab="", ylab="")
        rasterImage(img,0,0,1,1)
      }
      dev.off()
      a<-apply(as.data.frame(files),1,function(x) unlink(x))
    }  
  } 
}
project<-function(){
  path=paste(getwd(),"/",label,"/",sep="")
  system(paste("mkdir -p ", path,sep=""))
  cat(paste("All files will be saved under folder:", path),"\n")
  return(path)
}
