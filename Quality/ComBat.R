#!/usr/bin/env Rscript

library(limma)
library(sva)
library(edgeR)
library(usethis)
library(Glimma)

combat_analysis<-function(data,targetinfo,batch,group) {
  rawdata<-data
  colnames(rawdata)<-gsub("_nat","",colnames(rawdata))
  head(rawdata)
  targets<-targetinfo
  targets$Filename<-as.vector(targets$Filename)
  counts <- as.matrix(rawdata)
  batch <- targets$Patients
  group<-targets$Type
  counts_matrix <- as.matrix(counts)
  adjusted_counts <- ComBat_seq(counts, batch=batch, group=group, full_mod=FALSE)
  write.table(adjusted_counts,file="comBat_matrix.tab",quote = F,sep="\t",row.names =T)
}
