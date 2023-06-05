library(xlsx)
library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(ggpubr)
library(scales)

source("WCM_NEPC_surface_target_pipeline/gene_level_GTEx_helper_functions.R")
load("data/surface_protein_lists.RData")

#########
# Input #
#########

# benign.TPM : matrix (rows genes, columns subject) contain gene TPM of benign prostate samples
# PCa.TPM  : matrix (rows genes, columns subject) contain gene TPM of PCa samples
# CRPC.TPM  : matrix (rows genes, columns subject) contain gene TPM of CRPC samples
# NEPC.TPM  : matrix (rows genes, columns subject) contain gene TPM of NEPC samples

# GTEx.TPM  : matrix (rows genes, columns subject) contain gene TPM of GTEx tissue samples
# GTEx.benign.TPM  : matrix (rows genes, columns subject) contain gene TPM of GTEx normal prostate samples
# GTEx_sample_tissue : list corresponding to tissue types for each sample in GTEX.isoform
# GTEx_sample_tissue_subtype2 : list corresponding to tissue subtypes for each sample in GTEX.isoform

# NEPCvBenign : DESeq2 result data frame of NEPC compared to benign prostate samples
# NEPCvPCa : DESeq2 result data frame of NEPC compared to PCa samples
# NEPCvCRPC : DESeq2 result data frame of NEPC compared to CRPC samples

#redflag_tissues <- table containing tissue importance labels for each tissue in GTEx ( 0 = dispensible, 1 = non-critical, 2 = critical )

save_results <- FALSE # indicates whether to save result tables and plots

# Identify differentially expressed genes associated wth NEPC
NEPCvBenign_degs <- NEPCvBenign %>% filter(padj < 0.05 & log2FoldChange > log2(2)) %>% pull(Gene); length(NEPCvBenign_degs)
NEPCvPCa_degs <- NEPCvPCa %>%filter(padj < 0.1 & log2FoldChange > log2(1.5)) %>% pull(Gene); length(NEPCvPCa_degs)
NEPCvCRPC_degs <- NEPCvCRPC %>%filter(padj < 0.1 & log2FoldChange > log2(1.5)) %>% pull(Gene); length(NEPCvCRPC_degs)

#candidate_gene_list <- union(NEPCvBenign_degs, NEPCvPCa_degs)
candidate_gene_list <- intersect(NEPCvBenign_degs, NEPCvPCa_degs)
length(candidate_gene_list)

# Count the number of times that a gene is expressed in normal tissues of each importance type
gtex_expr_cutoff <- 3
candidate_GTEX_analysis <- candidates_GTEX_analysis(union(NEPCvBenign_degs, NEPCvPCa_degs))
GTEx.analysis2 <- candidates_GTEX_analysis_subtype(union(NEPCvBenign_degs, NEPCvPCa_degs))
candidate_GTEx_redflag_analysis <- data.table(
  Gene = rownames(GTEx.analysis2),
  GTEX.count = rowSums(candidate_GTEX_analysis > gtex_expr_cutoff),
  GTEX.count_critical = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==2,1])]>2),
  GTEX.count_noncritical = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==1,1])]>2),
  GTEX.count_disposable = rowSums(GTEx.analysis2[,which(colnames(GTEx.analysis2) %in% Janssen_redflag[Janssen_redflag[,2]==0,1])]>2)
)


# Form candidate list for putative NEPC-specific targets
candidate_df <- merge(
  NEPCvBenign %>% select(Gene, log2FoldChange, pvalue, padj),
  NEPCvPCa %>% select(Gene, log2FoldChange, pvalue, padj),
  by = "Gene",
  suffixes = c(".NvB", "")
) %>%
  merge(
    NEPCvCRPC %>% select(Gene, log2FoldChange, pvalue, padj),
    by = "Gene",
    suffixes = c(".NvP", ".NvC")
  ) %>%
  merge(
    candidate_GTEx_redflag_analysis,
    by = "Gene"
  ) %>%
  filter(Gene %in% union(NEPCvBenign_degs, NEPCvPCa_degs) & Gene %in% rownames(candidate_GTEX_analysis)) %>%
  mutate(
    medianBenign = sapply(Gene,function(this.gene) median(as.numeric(benign.TPM[rownames(benign.TPM)==this.gene,]))),
    medianPCa = sapply(Gene,function(this.gene) median(as.numeric(PCa.TPM[rownames(PCa.TPM)==this.gene,]))),
    medianCRPC = sapply(Gene,function(this.gene) median(as.numeric(CRPC.TPM[rownames(CRPC.TPM)==this.gene,]))),
    medianNEPC = sapply(Gene,function(this.gene) median(as.numeric(NEPC.TPM[rownames(NEPC.TPM)==this.gene,]))),
    surface = Gene %in% c(CD.surface.proteins, CSPA_hiconf,CSPA_hiconf1,CSPA_putative,CD.surface.proteins, GO0009897),
    surface_loose = Gene %in% c(CD.surface.proteins, CSPA_hiconf,CSPA_hiconf1,CSPA_putative,CD.surface.proteins, GO0009897, ion.goids, transmembrane.goids)
  )


general_candidate_list <- candidate_df %>% 
  filter(
    medianBenign < 1 &
      medianPCa < 1 & 
      medianNEPC > 1 & 
      surface_loose &
      GTEX.count <= 5
    ) %>% pull(Gene); length(general_candidate_list)
NEPC_candidate_list <- intersect(general_candidate_list, NEPCvCRPC_degs); length(NEPC_candidate_list)


if(save_results){
  candidate_df %>%
    filter(Gene %in% general_candidate_list) %>%
    mutate(
      `Candidate Type` = if_else(Gene %in% NEPC_candidate_list, "NEPC-specific", "NEPC+CRPC")
    ) %>%
    dplyr::rename_with(~ gsub("NvP", "NEPCvPCA", gsub("NvC", "NEPCvCRPC", gsub("NvB", "NEPCvBenign", .x, fixed = TRUE)))) %>%
    dplyr::rename(surface_strict = surface, surface = surface_loose) %>%
    write.xlsx("Top Candidate List.xlsx", sheetName = "Gene-Level")
}


# Heatmap of top candidates

mat <- candidate_df %>% 
  rename_with(~ gsub("GTEX.count_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("median", "", .x, fixed = TRUE)) %>%
  mutate(GTEx = critical/10) %>%
  select(Gene, GTEx, Benign, PCa, CRPC, NEPC) %>%
  filter(Gene %in% general_candidate_list) %>%
  tibble::column_to_rownames("Gene")

hm <- Heatmap(
  mat,
  name = "Expression\n# Tissues/\nTPM",
  row_names_gp = gpar(fontsize = 4, fontface = 2),
  row_split = factor(if_else(rownames(mat) %in% NEPC_candidate_list, "NEPC-specific\nCandidates", "NEPC+CRPC\nCandidates"), levels = c("NEPC+CRPC\nCandidates", "NEPC-specific\nCandidates")),
  col = circlize::colorRamp2(c(0,5), c("white", "navyblue")),
  cluster_columns = FALSE, 
  cluster_column_slices = FALSE,
  cluster_row_slices = FALSE,
  row_order = order(mat$CRPC, decreasing = TRUE),
  column_split = factor(c("GTEx", "Prostate", "PCA", "CRPC", "NEPC"), levels = c("GTEx", "Prostate", "PCA", "CRPC", "NEPC")),
  column_gap = unit(0, "mm"),
  border_gp = gpar(col = "black", lwd = 2.5),
  column_title_gp = gpar(fill = c("black", "white", "yellow", "forestgreen", "red"), fontface = 2, col = c("white", rep("black", 4)), fontsize = 10),
  column_title_side = "bottom",
  show_column_names = FALSE, 
  border = TRUE, 
  show_heatmap_legend = FALSE
); draw(hm, merge=T)

if(save_plots){
  pdf('Top Candidates Heatmap.pdf', width = 4, height = 7)
  draw(hm, merge=T)
  dev.off()
}



# Plot boxplots of individual gene expression patterns
gene_shortlist <- c("CEACAM5", "CELSR3")
tpm_df <- data.table()
for(gene in gene_shortlist){
  tpm_df <- rbind(
    tpm_df,
    data.table(
      Gene = gene,
      TPM = unlist(as.vector(GTEx.benign.TPM[which(rownames(GTEx.benign.TPM) == gene),])),
      type = "GTEx"
    ),
    data.table(
      Gene = gene,
      TPM = unlist(as.vector(benign.TPM[which(rownames(benign.FPKM) == gene),])),
      type = "Prost."
    ),
    data.table(
      Gene = gene,
      TPM = unlist(as.vector(PCa.TPM[which(rownames(PCa.TPM) == gene),])),
      type = "PCA"
    ),
    data.table(
      Gene = gene,
      TPM = unlist(as.vector(CRPC.TPM[which(rownames(CRPC.TPM) == gene),])),
      type = "CRPC"
    ),
    data.table(
      Gene = gene,
      TPM = unlist(as.vector(NEPC.TPM[which(rownames(NEPC.TPM) == gene),])),
      type = "NEPC"
    )
  )
}

p <- tpm_df %>%
  mutate(
    logTPM = log10(TPM + 1),
    gene = factor(gene, levels = gene_shortlist)
  ) %>%
  ggboxplot(
    x = "type",
    y = "logTPM",
    fill = "type",
    facet.by = "Gene",
    scales = "free_y",
    ncol = 1
  ) +
  scale_fill_manual(values = c("grey50", "white", "#CCCC00", "forestgreen", "#CB0202")) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = log10(1+c(0, 1, 10, 100, 1000)), labels = c(0, 1, 10, 100, 1000)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(
    c("NEPC", "CRPC"), 
    c("NEPC", "PCA"), 
    c("NEPC", "Prost."), 
    c("NEPC", "GTEx") ,
    c("CRPC", "PCA"), 
    c("CRPC", "Prost."), 
    c("CRPC", "GTEx") 
  )) +
  labs(
    x = "", 
    y = "TPM"
  )


if(save_plots){
  ggsave(
    filename = "Top Candidates Boxplots.pdf", 
    plot = p,
    width = 3.75,
    height = 7, 
    useDingbats = FALSE
  )
}

