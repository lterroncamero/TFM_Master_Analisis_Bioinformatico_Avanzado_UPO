#!/usr/bin/env Rscript

 #Load libraries
 library("limma")
 library("HsAgilentDesign026652.db")
 library("annotate")
 library("AnnotationDbi")
 library("org.Hs.eg.db")
 library ("ggrepel")

array_agilent_matrix<-function(targetinfo,label){
  #Set-up
  targetinfo <- readTargets("targetinfo",row.names="FileName",sep="")
  #Read files
  project <- read.maimages(targetinfo,source="genepix", green.only=TRUE)
  #Background correction
  project.bgc <- backgroundCorrect(project, method="normexp", offset=50)
  #Normalize the data with the 'quantile' method for 1-color
  project.NormData <-normalizeBetweenArrays(project.bgc,method="quantile")
  write.table(project.NormData$E, paste("Matriz_nor",label,".xls",sep=""), sep="\t", quote=FALSE, col.names = NA)
}

comparative_array<-function(data,treatment,level1,level2,covariate){
  #Create model
  treatment <- factor(treatment, levels=c(level1,level2))
  covariate<-factor(covariate)
  design <- model.matrix(~treatment+covariate)
  #adding gene names & colapsing
  df<-as.data.frame(project.NormData$E)
  df$Probes<-project.NormData$genes$ID
  df_agg<-aggregate(. ~ Probes, data = df, mean,na.rm=T)
  df_agg$Genes<-getSYMBOL(as.character(df_agg$Probes), "HsAgilentDesign026652")
  df_agg<-aggregate(df_agg[ , 2:11], by = list(df_agg$Genes), FUN = mean, na.rm=TRUE)
  rownames(df_agg)<-df_agg$Genes
  df_agg<-df_agg[setdiff(names(df_agg),"Genes")]
  fit <- lmFit(df_agg, design)
  fit2 <- eBayes(fit)
  topTable(fit2, coef="fLH7")
  probeset.list<-topTable(fit2, coef="fLH7", number=nrow(fit), adjust.method="BH", sort.by="P")
  results <- decideTests(fit2)
  summary(results)
  x<-as.data.frame(mapIds(
    org.Hs.eg.db, 
    keys = rownames(probeset.list), 
    "GENENAME", 
    "SYMBOL",
    fuzzy = TRUE,
    multiVals = "first"))
  probeset.list$Description<-x[,1]
  
  rownames(probeset.list) <- NULL
  colnames(probeset.list)[which(names(probeset.list) == "Group.1")] <- "Gene"
  write.table(probeset.list, paste("Matriz_nor",level1,"vs",level2,".xls",sep=""), sep="\t", quote=FALSE,col.names = NA)
}
  
Volcano_plot<-function(data,file,main){ 
 data$Gene<-rownames(data)
 data$Threshold<-"nDEG"
 data[data$adj.P.Val <= 0.05,"Threshold"]<-"DEG"
 ##Top20
 real_DE<-data[data$adj.P.Val<=0.05,]
 selected_FC<-head(real_DE[order(abs(real_DE$logFC),decreasing = T),], 20)
 mycolors <- c("purple", "gray")
 names(mycolors) <- c("DEG", "nDEG")
 pdf("Vocano_plot_LH7_vs_LH2_purple.pdf",paper="a4")
  g = ggplot(data=data, aes(x=logFC, y=-log10(adj.P.Val), colour=Threshold)) + coord_cartesian(xlim = c(-6, 6 )) +
  geom_point(alpha=0.4, size=1.75) +
  xlab("log2 fold change") + ylab("-log10 adj.P.Val") + theme_minimal() + scale_colour_manual(values = mycolors) +
  geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3) + ggtitle("LH+7 vs LH+2") + theme(plot.title = element_text(hjust = 0.5)) #adding text for the top1 20 genes
  plot(g)
 dev.off()
 png("Vocano_plot_LH7_vs_LH2_purple.png", res=300, width=2000, height=2000)
  g = ggplot(data=data, aes(x=logFC, y=-log10(adj.P.Val), colour=Threshold)) + coord_cartesian(xlim = c(-6, 6 )) +
  geom_point(alpha=0.4, size=1.75) +
  xlab("log2 fold change") + ylab("-log10 P.Val") + theme_minimal() + scale_colour_manual(values = mycolors) +
  geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3) + ggtitle("LH+7 vs LH+2") + theme(plot.title = element_text(hjust = 0.5)) #adding text for the top1 20 genes
  plot(g)
 dev.off()
}


