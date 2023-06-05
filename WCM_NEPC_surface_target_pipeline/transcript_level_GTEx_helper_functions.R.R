candidates_GTEX_analysis <- function(signature){
  GTEx.rpkm.surface <- GTEX.isoform[gsub("[.].*$","",rownames(GTEX.isoform)) %in% signature, order(GTEx_sample_tissue)]
  
  if(is.null(dim(GTEx.rpkm.surface))){
    GTEx.rpkm.surface.2 <- matrix(GTEx.rpkm.surface[order(names(GTEx.rpkm.surface))],nrow=1)
  }else{
    GTEx.rpkm.surface.2 <- GTEx.rpkm.surface[order(rownames(GTEx.rpkm.surface)),]
  }
  
  GTEx.analysis <- c()
  for(tissue in unique(GTEx_sample_tissue)){
    temp <- GTEx.rpkm.surface.tissue <- GTEx.rpkm.surface.2[,GTEx_sample_tissue2==tissue];
    GTEx.rpkm.surface.tissue <- matrix(as.numeric(unlist(GTEx.rpkm.surface.tissue)), dim(temp));
    rownames(GTEx.rpkm.surface.tissue) <- rownames(temp);
    colnames(GTEx.rpkm.surface.tissue) <- colnames(temp)
    
    temp <- apply(GTEx.rpkm.surface.tissue,1,median)
    names(temp) <- rownames(GTEx.rpkm.surface.tissue)
    
    if(ncol(GTEx.rpkm.surface.tissue) > 0){
      GTEx.analysis <- cbind(GTEx.analysis,temp)
      colnames(GTEx.analysis)[ncol(GTEx.analysis)] <- tissue    
    }
  }
  
  return(GTEx.analysis)
}

candidates_GTEX_analysis_subtype <- function(total.sig.surface){
  GTEx.rpkm.surface <- GTEX.isoform[gsub("[.].*$","",rownames(GTEX.isoform)) %in% total.sig.surface, order(GTEx_sample_tissue_subtype)]
  GTEx_sample_tissue_subtype2 <- GTEx_sample_tissue_subtype[order(GTEx_sample_tissue_subtype)]
  
  if(is.null(dim(GTEx.rpkm.surface))){
    GTEx.rpkm.surface.2 <- matrix(GTEx.rpkm.surface[order(names(GTEx.rpkm.surface))],nrow=1)
  }else{
    GTEx.rpkm.surface.2 <- GTEx.rpkm.surface[order(rownames(GTEx.rpkm.surface)),]
  }
  
  GTEx.analysis <- c()
  for(subtissue in unique(GTEx_sample_tissue_subtype)){
    GTEx.rpkm.surface.tissue <- GTEx.rpkm.surface.2[,GTEx_sample_tissue_subtype2==subtissue];
    temp <- GTEx.rpkm.surface.tissue; 
    GTEx.rpkm.surface.tissue <- matrix(as.numeric(unlist(GTEx.rpkm.surface.tissue)),dim(temp))
    rownames(GTEx.rpkm.surface.tissue) <- rownames(temp);
    colnames(GTEx.rpkm.surface.tissue) <- colnames(temp)
    
    temp <- as.numeric(apply(GTEx.rpkm.surface.tissue,1,median,na.rm=T))
    names(temp) <- rownames(GTEx.rpkm.surface.tissue)
    
    if(ncol(GTEx.rpkm.surface.tissue)  > 0){
      GTEx.analysis <- cbind(GTEx.analysis, temp)
      colnames(GTEx.analysis)[ncol(GTEx.analysis)] <- subtissue    
    }
  }
  colnames(GTEx.analysis) <- gsub("[ ]\\(.*$","",(colnames(GTEx.analysis)))

  return(GTEx.analysis)
}