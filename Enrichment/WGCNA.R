#!/usr/bin/env Rscript

library("WGCNA")

WGCNA<-function(nexpr_mat, minSize=30, MEDissThres=0.25, organism) {
  
  nexpr_mat<-df_agg
  minSize<-30
  MEDissThres<-0.40
  organism<-"human"
  
  # min(20, ncol(datExpr)/2 )
  # library(WGCNA)
  # library(tidyr)
  
  allowWGCNAThreads()
  
  path<-"/mnt/beegfs/Bioinformatica/EXT/08-EXT-20-Saltmae-endometrio/Analisis/arrays/"
  new_path=paste(path,"WGCNA_sin_corregir/",sep="")
  system(paste("mkdir -p ", new_path,sep=""))
  
  #rownames(nexpr_mat)<-nexpr_mat$Group.1
  #nexpr_mat <- nexpr_mat[,-1]
  
  t_exprs=t(nexpr_mat)
  
  ## Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(t_exprs, powerVector = powers)
  if (is.na(sft$powerEstimate)){
    power = sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
  }  else {power = sft$powerEstimate}
  
  pindex=which(sft$fitIndices$Power==power)
  
  # Plot results
  label<-"_no_corregidos_0.50"
  pdf(file =  paste(new_path,"WGCNA_plots_",label,".pdf",sep=""),paper="a4")
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[-pindex,1], -sign(sft$fitIndices[-pindex,3])*sft$fitIndices[-pindex,2],
       labels=powers[-pindex],col="black", cex=0.7)
  text(sft$fitIndices[pindex,1], -sign(sft$fitIndices[pindex,3])*sft$fitIndices[pindex,2],
       labels=power,col="red")
  abline(h=0.90,col="blue")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[-pindex,1], sft$fitIndices[-pindex,5], labels=powers[-pindex],col="black")
  text(sft$fitIndices[pindex,1], sft$fitIndices[pindex,5], labels=power,col="red")
  par(mfrow = c(1,1))
  
  ## Correlation network construction and module identification (one-step)
  temp_cor <- cor       
  cor <- WGCNA::cor
  
  sink(file = paste0(new_path,label,"_blockwiseModules.out"))
  netwk <- blockwiseModules(t_exprs,
                            corType = "pearson",
                            
                            # == Adjacency Function ==
                            power = power,
                            networkType = "signed",
                            
                            # == Tree and Block Options ==
                            deepSplit = 2,              # 0 - 4. simplified control over how sensitive module detection should be to module splitting
                            pamRespectsDendro = F,
                            # detectCutHeight = 0.75,
                            minModuleSize = minSize,
                            
                            # == Module Adjustments ==
                            reassignThreshold = 0,
                            mergeCutHeight = MEDissThres,
                            
                            # == TOM == Archive the run results in TOM file (saves time)
                            # saveTOMs = T,
                            # saveTOMFileBase = paste0(label,"_blockwiseTOM"),
                            
                            # == Output Options
                            numericLabels = F,
                            verbose = 3)
  sink()
  
  cor <- temp_cor
  
  ## Plot module merge results
  # Calculate dissimilarity of module eigengenes
  MEList = moduleEigengenes(t_exprs, colors = netwk$unmergedColors)
  MEs = MEList$eigengenes
  colnames(MEs)=gsub("ME","",colnames(MEs))
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  plot(METree, main = paste0("Module eigengenes merge result (",
                             length(unique(netwk$unmergedColors))," -> ",
                             length(unique(netwk$colors)),")"),
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  
  
  # Plot the dendrogram and the module colors underneath
  for (i in 1:length(netwk$dendrograms)) {
    plotDendroAndColors(netwk$dendrograms[[i]],
                        netwk$colors[netwk$blockGenes[[i]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05 )
  }
  
  # Get Module Eigengenes per cluster
  MEs0 = orderMEs(netwk$MEs)
  # Reorder modules so similar modules are next to each other
  module_order = gsub("ME","",names(MEs0))
  
  # Add treatment names
  MEs0$treatment = group
  
  # tidy & plot data
  mME = MEs0 %>%
    pivot_longer(-treatment) %>%
    mutate(name = gsub("ME", "", name), name = factor(name, levels = module_order))
  
  mt.plot <-ggplot(mME, aes(x=treatment, y=name, fill=value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1,1)) +
    theme(axis.text.x = element_text(angle=90, size = 8, vjust = 0.5, hjust = 1)) +
    labs(title = "Module-type Relationships", y = "Modules", x="Condition", fill="corr")
  plot(mt.plot)
  dev.off()
  
  jpeg(paste(new_path,"WGCNA_module_type_",label,".jpeg",sep=""))
  plot(mt.plot)
  dev.off()
  
  # Relate Module to Type Groups
  module_df = data.frame(
    gene_id = names(netwk$colors),
    colors = netwk$colors)
  write.table(module_df, file = paste(new_path,"WGCNA_modules_",label,".tsv",sep=""),sep="\t",row.names = F)
  
  group_enrichment(genes=names(netwk$colors) , mergedColors=netwk$colors, organism= organism, new_path=new_path)
  
}
group_enrichment <- function(genes, mergedColors, organism, ontology="BP, CC, MF",new_path){
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
  
  entrezid = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=dbName)
  
  GOenr = GOenrichmentAnalysis(mergedColors, entrezid$ENTREZID, organism = organism, nBestP = 10, getTermDetails = F,
                               verbose = 0) #entrezid, organism
  
  tab = GOenr$bestPTerms[[ontology]]$enrichment[, c(4,1,11, 5, 6, 12, 13, 2, 7)] #ontology ("BP" "CC" "MF" o "BP, CC, MF")
  
  # Round the numeric columns to 2 decimal places:
  tab[, c(4, 5)] = signif(apply(tab[, c(4, 5)], 2, as.numeric), 2)
  
  # Add genenames to GOTerms
  names(GOenr$bestPTerms[[ontology]]$forModule)=unique(tab$module)
  # tab$gene_names=apply(tab,1,FUN = function(x){ paste(entrezid[GOenr$bestPTerms[[ontology]]$forModule[[x[2]]][[as.integer(x[1])]]$genePositions,"SYMBOL"], collapse = ", ") })
  
  tab=tab[,-1]
  colnames(tab) = c("module","termID", "p-val", "adj-pval", "ont", "term name", "modSize", "nInTerm") # , "Gene_names"
  write.table(tab, file = paste(new_path,"WGCNA_GOEnrichment_",label,".tsv",sep=""), sep = "\t", quote = TRUE, row.names = FALSE)
}
plot_module <- function(colour, nexpr_mat){
  
  colour<-"blue"
  nexpr_mat<-df_agg
  
  # Pull out list of genes in that module
  submod = read.delim(paste0(path,"WGCNA/WGCNA_modules_",label,".tsv")) %>% subset(colors == colour)
  
  # Get normalized expression for those genes
  subexpr = nexpr_mat[submod$gene_id,]
  colnames(subexpr)= targets$FileName
  
  submod_df = data.frame(subexpr) %>%
    mutate(
      gene_id = row.names(.)
    ) %>%
    pivot_longer(-gene_id)
  
  submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
    geom_line(aes(alpha = 0.2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 ,vjust = 0.5, size=10), legend.position = "None") +
    labs(x = "Sample",
         y = "Normalized expression",
         title = paste("Expression of genes from",colour,"module across samples"))

}
