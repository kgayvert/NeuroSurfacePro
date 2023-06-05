library(data.table)
library(gplots)
library(ggplot2)
library(plyr)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(GenomicRanges)
library(ggpubr)
library(scales)

source("WCM_NEPC_surface_target_pipeline/transcript_level_GTEx_helper_functions.R")
load("data/surface_protein_lists.RData")

#########
# Input #
#########

# benign.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of benign prostate samples
# PCa.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of PCa samples
# CRPC.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of CRPC samples
# NEPC.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of NEPC samples
# H660.FPKM :  matrix (rows genes, columns subject) contain isoform FPKM of H660 cell line

# GTEX.isoform : matrix (rows genes, columns subject) contain isoform FPKM of GTEx tissue samples
# GTEx_sample_tissue : list corresponding to tissue types for each sample in GTEX.isoform
# GTEx_sample_tissue_subtype2 : list corresponding to tissue subtypes for each sample in GTEX.isoform

# NEPCvBenign : isoform-level cuffdiff result data frame of NEPC compared to benign prostate samples
# NEPCvPCa : isoform-level cuffdiff result data frame of NEPC compared to PCa samples
# NEPCvCRPC : isoform-level cuffdiff result data frame of NEPC compared to CRPC samples

#redflag_tissues <- table containing tissue importance labels for each tissue in GTEx ( 0 = dispensible, 1 = non-critical, 2 = critical )

save_results <- FALSE # indicates whether to save result tables and plots

load("../data/surface_protein_lists.RData")

# Add annotation on whether transcript is thought to be protein coding
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
protein_coding_status_gb <- getBM(attributes=c("ensembl_transcript_id","transcript_biotype"),filters = c("ensembl_transcript_id"), values=list(NEPCvBenign$Isoform), mart=ensembl)
protein_coding_transcripts <- protein_coding_status_gb %>% filter(transcript_biotype == "protein_coding") %>% pull(ensembl_transcript_id) 

NEPCvBenign <- NEPCvBenign %>%
  mutate(
    protein_coding= Isoform %in% protein_coding_transcripts 
)

# Start pipeline
fcthresh <- 2
pthresh <- 0.05
NEPCvBenign.sig <- NEPCvBenign[NEPCvBenign$p_value < pthresh & NEPCvBenign$log2.fold_change. > log2(fcthresh),] %>% 
  mutate(
    surface = gene %in% c(CD.surface.proteins, CSPA_hiconf,CSPA_hiconf1,CSPA_putative,CD.surface.proteins, GO0009897),
    surface_loose = gene %in% c(CD.surface.proteins, CSPA_hiconf,CSPA_hiconf1,CSPA_putative,CD.surface.proteins, GO0009897, ion.goids, transmembrane.goids)
  ) %>%
  dplyr::select(Isoform, gene, log2.fold_change., p_value, protein_coding, surface, surface_loose)

# Count the number of times that a gene is expressed in normal tissues of each importance type
gtex_expr_cutoff <- 3
candidate_GTEX_analysis <- candidates_GTEX_analysis_v2(NEPCvBenign.sig$Isoform)
GTEx.analysis2 <- candidates_GTEX_analysis_subtype(NEPCvBenign.sig$Isoform,"all",plot=FALSE)$GTEx.analysis

# Form candidate list for putative NEPC-specific targets
candidate_isoform_df <- NEPCvBenign.sig %>%
  filter(Isoform %in% rownames(GTEx.analysis2)) %>%
  mutate(
    GTEX.count = rowSums(candidate_GTEX_analysis > gtex_expr_cutoff),
    GTEX.count_critical = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==2,1])]>2),
    GTEX.count_noncritical = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==1,1])]>2),
    GTEX.count_disposable = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==0,1])]>2),
    medianNEPC = sapply(Isoform,function(x){val=NA;try(silent = TRUE,expr={val=median(as.numeric(NEPC.iso.FPKM[grep(x,NEPC.iso.FPKM$Transcript),3:ncol(NEPC.iso.FPKM)]),na.rm=T)});return(val)}),
    medianCRPC = sapply(Isoform,function(x){val=NA;try(silent = TRUE,expr={val=median(as.numeric(CRPC.iso.FPKM[grep(x,CRPC.iso.FPKM$Transcript),3:ncol(CRPC.iso.FPKM)]),na.rm=T)});return(val)}),
    medianPCa = sapply(Isoform,function(x){val=NA;try(silent = TRUE,expr={val=median(as.numeric(PCa.WCMC.iso.FPKM[grep(x,PCa.WCMC.iso.FPKM$Transcript),3:ncol(PCa.WCMC.iso.FPKM)]),na.rm=T)});return(val)}),
    medianBenign = sapply(Isoform,function(x){val=NA;try(silent = TRUE,expr={val=median(as.numeric(benign.WCMC.iso.FPKM[grep(x,benign.WCMC.iso.FPKM$Transcript),3:ncol(benign.WCMC.iso.FPKM)]),na.rm=T)});return(val)})
  )

# Generate candidate list
candidate_list <- candidate_isoform_df %>% 
  filter(
    !is.na(protein_coding) & protein_coding &
      GTEX.count == 0 &
      medianBenign < 1 &
      medianPCa < 1 & 
      medianNEPC > 1 & 
      surface_loose
  ) %>%
  pull(Isoform)

NEPC_candidate_list <- candidate_isoform_df %>% 
  filter(
    GTEX.count == 0 &
      medianBenign < 1 &
      medianPCa < 1 & 
      medianCRPC < 1 & 
      medianNEPC > 2 & 
      medianNEPC / medianCRPC > log2(1.5) &
      surface_loose
  ) %>% 
  pull(Isoform)

if(save_results){
  candidate_isoform_df %>%
    filter(Isoform %in% candidate_list) %>%
    mutate(
      `Candidate Type` = if_else(Isoform %in% NEPC_candidate_list, "NEPC-specific", "NEPC+CRPC")
    ) %>%
    dplyr::rename(
      Transcript = Isoform, 
      Gene = gene,
      log2FoldChange.NEPCvBenign = log2.fold_change.,
      pvalue.NEPCvBenign = p_value,
      surface_strict = surface, 
      surface = surface_loose
    ) %>%
    dplyr::select(-protein_coding, -surface) %>%
    dplyr::select(!matches("GTEx")) %>%
    write.xlsx("Top Candidate List.xlsx", sheetName = "Transcript-Level", append = TRUE, row.names = FALSE)
  
}

# Plot boxplots of individual transcript expression patterns
select_transcripts <- c("ENST00000399380", "ENST00000269141")
tpm_df <- data.table()
for(isoform in select_transcripts){
  tpm_df <- rbind(
    tpm_df,
    data.table(
      Isoform = isoform,
      FPKM = as.numeric(unlist(as.vector(GTEX.isoform[which(rownames(GTEX.isoform) == isoform),grep("Prostate", count.sample.tissue)]))),
      type = "GTEx"
    ),
    data.table(
      Isoform = isoform,
      FPKM = as.numeric(unlist(as.vector(benign.WCMC.iso.FPKM[which(benign.WCMC.iso.FPKM[,1] == isoform),-1*1:2]))),
      type = "Prost."
    ),
    data.table(
      Isoform = isoform,
      FPKM = as.numeric(unlist(as.vector(PCa.WCMC.iso.FPKM[which(PCa.WCMC.iso.FPKM[,1] == isoform),-1*1:2]))),
      type = "PCA"
    ),
    data.table(
      Isoform = isoform,
      FPKM = as.numeric(unlist(as.vector(CRPC.iso.FPKM[which(CRPC.iso.FPKM[,1] == isoform),-1*1:2]))),
      type = "CRPC"
    ),
    data.table(
      Isoform = isoform,
      FPKM = as.numeric(unlist(as.vector(NEPC.iso.FPKM[which(NEPC.iso.FPKM[,1] == isoform),-1*1:2]))),
      type = "NEPC"
    )
  )
}

p <- tpm_df %>%
  mutate(
    logFPKM = log10(FPKM + 1),
    Isoform = factor(Isoform, levels = select_transcripts)
  ) %>%
  ggboxplot(
    x = "type",
    y = "logFPKM",
    fill = "type",
    facet.by = "Isoform",
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(values = c("grey50", "white", "#CCCC00", "forestgreen", "#CB0202")) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = log10(1+c(0, 1, 10, 100, 1000)), labels = c(0, 1, 10, 100, 1000)) +
  labs(
    x = "", 
    y = "FPKM"
  )

if(save_plots){
  ggsave(
    filename = "Isoform Boxplots.pdf", 
    plot = p, 
    width = 7, 
    height = 3,
    useDingbats = FALSE
  )
}
