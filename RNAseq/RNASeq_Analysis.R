
library("edgeR")
library("RColorBrewer")
library("genefilter")
library("ggrepel")
library("ggfortify")
library("cluster")
library("factoextra")
library("ggplot2")
library("M3C")
library("AnnotationDbi")
library("xlsx")
library("clusterProfiler")
library("png")
library("curl")
library("corrplot")
library("tidyr")
library("dplyr")
library("WGCNA")

Lau_QC<-function(label,dge,data){
  x<-dge
  group<-x$samples$group
  
  log_Edu<-log(rawdata,2)
  cpm <- cpm(x)
  lcpm <- cpm(x, log=TRUE)
  
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  c(L, M)
  
  summary(lcpm)
  table(rowSums(x$counts==0)==6)
  
  pdf(paste(path,label,"_QC.pdf",sep=""),paper="A4")
  lcpm.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(log_Edu[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(log_Edu[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright",legend=targets$Name, text.col=col, bty = "n", cex = 0.5)
  
  lcpm <- cpm(x, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=targets$Name, text.col=col, bty="n", cex = 0.5)
  
  x <- calcNormFactors(x)
  x$samples
  x <- estimateCommonDisp(x, robust=TRUE)
  x <- estimateTagwiseDisp(x)
  
  x2 <- x
  x2$samples$norm.factors <- 1

  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  
  par(mfrow=c(1,2))
  lcpm <- cpm(dge, log=TRUE)
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Unnormalised data",ylab="Log-cpm")
  
  lcpm <- cpm(x, log=TRUE)
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Normalised data",ylab="Log-cpm")
  
  par(mfrow=c(1,1))
  
  #Libray size
  barplot(x$samples$lib.size,names.arg = targets$Name,las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  
  #Corrplot

  lcpm <- cpm(x, log= FALSE)
  corrplot(cor(lcpm,method="spearman"), method='number',type = 'upper')
  corrplot(cor(lcpm,method="spearman"), order='AOE')
  
  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  levels(col.group)
  z<-plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F)
  edge<-sd(z$x)
  plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="A. MDS-PCoA Sample Names")
  
  lcpm <- cpm(x, log=TRUE)
  par(mfrow=c(1,1))
  #col.group <- as.factor(targets$Type)
  #levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  #col.group <- as.character(col.group)
  #levels(col.group)
  z<-plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F)
  edge<-sd(z$x)
  #cat(edge)
  plotMDS(lcpm, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="A. MDS-PCoA log2 Sample Names")
  
  
  #PCA por tipos
  data_pca<-as.matrix(x)
  data_pca<-as.data.frame(t(data_pca))
  rownames(data_pca)<-targets$Filename
  data_pca.PC = prcomp(data_pca)
  data_pca$Name<-targets$Name
  data_pca$Filename<-targets$Filename
  data_pca$Name<-targets$Name
  data_pca$Sex<-targets$Sex
  data_pca$Age<-targets$Age
  data_pca$VAS_Group<-targets$VAS_Group
  data_pca$TypeII<-targets$TypeII
  
  
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Name',xlim = c(-10,10)))
  
  rsd <- rowSds(as.matrix(x))
  sel <- order(rsd, decreasing=TRUE)[1:250]
  samplenames<-as.character(targets$Name)
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities",cexRow=0.01,cexCol=0.5,labCol=samplenames)
  
  #cluster
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  #pr.hc.c<- hclust(na.omit(dist(t(data))))
  #plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of ", label, sep=""), labels=targets$Filemane, cex=0.5)
  #Normalized clustering analysis plot
  #pr.hc.c<- hclust(na.omit(dist(t(dgenorm$counts))))
  pr.hc.c<- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of Normalized samples of ", label, sep=""), labels=targets$Name, cex=0.5)
  
  #tSNE
  #a<-tsne(x$counts,seed=100,labels=as.factor(targets$Type), perplex=perplex, legendtitle="Types",text=targets$Type ,dotsize=3, legendtextsize = 8) + ggtitle("Tsne") + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  #plot(a)  
  dev.off()
  
  jpeg(paste0(path,"PCoA_QC.jpg"))
  plotMDS(lcpm, labels=targets$Type, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  dev.off()
  
  jpeg(paste0(path,"Samples_HTree_QC.jpg"))
  plot(pr.hc.c, xlab="Sample Distance", labels=targets$Name, cex=0.5)
  dev.off()
  
}
lauVolcano<-function(data,file,main){
  options(ggrepel.max.overlaps = Inf)
  data<-data$table
  data$Gene<-rownames(data)
  data$threshold<-"nDEG"
  data[data$FDR <= 0.05,"threshold"]<-"DEG"
  #
  data[data$FDR==0,"adj.P.Val"]<-1e-318
  #
  
  ##Los 10 con mejor FC
  real_DE<-data[data$FDR<=0.05,]
  selected_FC<-head(real_DE[order(abs(real_DE$logFC),decreasing = T),], 20)
  max_value=max(abs(real_DE$logFC))
  pdf(paste(path,file,sep=""),paper="a4")
  g = ggplot(data=data, aes(x=logFC, y=-log10(FDR), colour=threshold)) + coord_cartesian(xlim = c(-max_value, max_value )) +
    geom_point(alpha=0.4, size=1.75) +
    xlab("log2 fold change") + ylab("-log10 adj.P.Val") +
    geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3) + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) #adding text for the top1 20 genes
  plot(g)
  dev.off()
  jpeg(paste0(path, gsub(".pdf",".jpg", file)))
  plot(g)
  dev.off()
}
createExcel<-function(data,file,organism="human"){
  dbName<-NULL
  if(organism == "human"){
    library("org.Hs.eg.db")
    dbName<-org.Hs.eg.db
  }else if (organism == "mouse"){
    library("org.Mm.eg.db")
    dbName<-org.Mm.eg.db
    
  }else{
    stop("At the moment only human and mouse are supported")
  }
  options(java.parameters = "-Xmx2048m")

  wb<-createWorkbook(type="xlsx")
  
  CellStyle(wb, dataFormat=NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
  
  TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
    Alignment(wrapText=FALSE, horizontal="ALIGN_CENTER") +
    Border(color="black", position=c("TOP", "BOTTOM", "LEFT", "RIGHT"), 
           pen=c("BORDER_THIN", "BORDER_THIN", "BORDER_THIN", "BORDER_THIN"))
  
  TextFormat = CellStyle(wb, dataFormat=DataFormat("@"))
  FormatList = list('1'=TextFormat, '2'=TextFormat,'3' = TextFormat,'4' = TextFormat)
  
  sheet1 <- xlsx::createSheet(wb, sheetName='DE Genes')
  
  x<-as.data.frame(mapIds(
    dbName, 
    keys = data$Gene, 
    "GENENAME", 
    "SYMBOL",
    fuzzy = TRUE,
    multiVals = "first"))
  
  data$Description<-x[,1]
  
  data<-data[,c(1,7,2:6)]
  xlsx::addDataFrame(data, sheet1, col.names=TRUE, colStyle=FormatList, 
                     row.names=FALSE, colnamesStyle = TABLE_COLNAMES_STYLE)
  
  addAutoFilter(sheet1, paste(LETTERS[1],LETTERS[ncol(data)],sep="-"))
  for (column in 1:ncol(data)) {
    #print(paste("LA columna ", column, " tiene un total de ",  max(nchar(na.omit(data[,column])))))
    size=max(nchar(na.omit(data[,column])))
    if(size <10) size<-15
    setColumnWidth(sheet1, column, size)
  }
  
  saveWorkbook(wb, file = paste(path,file,sep=""))
  write.table(data, paste0(path, gsub("DEG_","",gsub(".xlsx",".tsv",file))), sep = "\t", row.names = F)
}

filter<-function(filter="standard",data,min_group=3){
  if(filter == "standard"){
    keep<-rowSums(cpm(data)>1) >= min_group
    data<-data[keep,]
    data$samples$lib.size <-colSums(data$counts)
  }
  else if(filter == "bin"){
    keep<-rowSums(data$counts)>0
    data<-data[keep,]
    data$samples$lib.size <-colSums(data$counts)
  }
  else{
    stop("At the moment only bin/standard are supported")
  }
  return(data)
}
marker_exp<-function(GeneName, save=F){
  
  group<-dgenorm$samples$group
  col.group <- as.factor(group)
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  
  shit<-as.vector(dgenorm[GeneName,]$counts)
  raw<-barplot(shit,names.arg = targets$Name,las=2,col=col.group,main=paste(GeneName, " raw_count expression", sep=""),ylim=range(pretty(c(0, shit))))
  plot_rpkm<-barplot(rpkm[GeneName,],names.arg = targets$Name,las=2,col=col.group,main=paste(GeneName, " RPKM expression", sep=""),ylim=range(pretty(c(0, rpkm[GeneName,]))))
  plot(raw)
  plot(plot_rpkm)
  
  if (save) {
    pdf(file = paste(path,GeneName,"_Expression.pdf"))
    plot(raw)
    plot(plot_rpkm)
    dev.off()
    
    png(file = paste(path,GeneName,"_Expression.png"))
    plot(plot_rpkm)
    dev.off()
  }
}
myDGE<-function(filter,min_group=3){
  #creating DGE object
  
  dge<-DGEList(counts=rawdata, group=group, genes=results_counts) 
  dge$sample
  
  #filtering low expressed genes
  table(dge$samples$group)
  cat(paste("Number of genes before filter [ ",filter," ]: ", nrow(dge),sep=""),"\n")
  dge<-filter(filter=filter,dge,min_group = min_group)
  cat(paste("Number of genes after filter [ ",filter," ]: ", nrow(dge),sep=""),"\n")
  nrow(dge)
  
  #normalizing
  dgenorm <- calcNormFactors(dge)
  dgenorm <- estimateCommonDisp(dgenorm, robust=TRUE)
  dgenorm <- estimateTagwiseDisp(dgenorm)
  
  #Quality
  Lau_QC(label=label,dge=dge,data=rawdata)
  
  return(dgenorm)
}
norm_exp<-function(type="RPKM"){
  #saving RPKM
  rpkm<-rpkm(dgenorm,normalized.lib.sizes=TRUE)
  colnames(rpkm)<-rownames(dgenorm$samples)
  resultsfile<-file.path(paste(path,label,"_RPKM.xls", sep=""))
  write.table(rpkm, file=resultsfile, sep = "\t", col.names = NA , qmethod = "double") 
  return(rpkm)
}
DGE<-function(comp,organism, myFDR=0.05, myFC=0){
  et <- exactTest(dgenorm,pair = comp)
  #Extracting the statistical data order by p-value
  top1<-topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue") 
  
  summary(decideTests(et))
  nrow(top1$table[top1$table$FDR<=myFDR & abs(top1$table$logFC)>= myFC,])
  print(head(top1$table,10))
  
  myLabel1=paste(comp, collapse = '_vs_')
  myLabel2=paste(myLabel1, "FDR", myFDR, "FC", myFC, sep="_")
  createExcel(top1$table,paste("DEG_",myLabel1,".xlsx", sep=""),organism = organism)
  
  #my_enrichment(top1$table,FA_label=myLabel2,cutoff=0.05,organism = organism, FDR=myFDR,FC=myFC)
  lauVolcano(top1,paste("Volcano_plot_", myLabel1 , ".pdf", sep=""),myLabel1)
  return(top1)
}
project<-function(){
  path=paste(getwd(),"/",label,"/",sep="")
  system(paste("mkdir -p ", path,sep=""))
  cat(paste("All files will be saved under folder:", path),"\n")
  return(path)
}
