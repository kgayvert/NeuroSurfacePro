candidates_GTEX_analysis <- function(genelist){
  GTEx.TPM.surface <- GTEx.TPM[rownames(GTEx.TPM ) %in% genelist, order(GTEx_sample_tissue)]
  if(is.null(dim(GTEx.TPM.surface))){
    GTEx.TPM.surface <- matrix(GTEx.TPM.surface[order(names(GTEx.TPM.surface))],nrow=1)
  }else{
    GTEx.TPM.surface <- GTEx.TPM.surface[order(rownames(GTEx.TPM.surface)),]
  }
  
  GTEx.analysis <- c()
  for(tissue in setdiff(GTEx.tissues, "Prostate")){
    GTEx.TPM.surface.tissue <- GTEx.TPM.surface[,GTEx_sample_tissue2==tissue];
    temp <- apply(GTEx.TPM.surface.tissue,1,median)
    if(ncol(GTEx.TPM.surface.tissue) > 0){
      GTEx.analysis<- cbind(GTEx.analysis,temp)
      colnames(GTEx.analysis)[ncol(GTEx.analysis)] <- tissue    
    }
  }
  return(GTEx.analysis)
}

candidates_GTEX_analysis_subtype <- function(total.sig.surface){
  GTEx.rpkm.surface <- GTEx.rpkm.values[rownames(GTEx.rpkm.values) %in% total.sig.surface,order(GTEx_sample_tissue_subtype)]
  GTEx_sample_tissue_subtype2 <- GTEx_sample_tissue_subtype[order(GTEx_sample_tissue_subtype)]
  
  if(is.null(dim(GTEx.rpkm.surface))){
    GTEx.rpkm.surface.2 <- matrix(GTEx.rpkm.surface[order(names(GTEx.rpkm.surface))],nrow=1)
  }else{
    GTEx.rpkm.surface.2 <- GTEx.rpkm.surface[order(rownames(GTEx.rpkm.surface)),]
  }
  
  GTEx.analysis <- c()
  for(subtissue in unique(GTEx_sample_tissue_subtype)){
    GTEx.rpkm.surface.tissue=GTEx.rpkm.surface.2[,GTEx_sample_tissue_subtype2==subtissue];
    temp=apply(GTEx.rpkm.surface.tissue,1,median)
    
    if(ncol(GTEx.rpkm.surface.tissue) > 0){
      GTEx.analysis=cbind(GTEx.analysis,temp)
      colnames(GTEx.analysis)[ncol(GTEx.analysis)] <- subtissue    
    }
  }
  colnames(GTEx.analysis) <- gsub("[ ]\\(.*$","",(colnames(GTEx.analysis)))

  return(GTEx.analysis)
}
