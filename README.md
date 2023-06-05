# NeuroSurfacePro

Multi-step bioinformatics pipeline to identify NEPC-specific, over-expressed gene transcripts that encode surface proteins.

This work is further described in: [citation to be added upon publication]

## Data availability
RNAseq data was obtained from  a previously published (Beltran et al., 2016; Berger et al, 2019) clinical cohort of tissue from NEPC tumors (n = 27), benign prostate tissue (n = 66), locally advanced prostate cancer (PCA, n = 68), and mCRPC (n = 72).

The GTEx data used for the analyses described were obtained from the [GTEx Portal](https://gtexportal.org/home/) on 03/01/2015 and/or dbGaP accession number [phs000424.v5.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v5.p1) on 03/01/2015


## Pipeline expected input:

NEPCvBenign : gene-level DESeq2 or isoform-level cuffdiff result data frame of NEPC compared to benign prostate samples

NEPCvPCa : gene-level DESeq2 or isoform-level cuffdiff result data frame of NEPC compared to PCa samples

NEPCvCRPC : gene-level DESeq2 or isoform-level cuffdiff result data frame of NEPC compared to CRPC samples4


redflag_tissues : table containing tissue importance labels for each tissue in GTEx ( 0 = dispensible, 1 = non-critical, 2 = critical )


GTEx_sample_tissue : list corresponding to tissue types for each sample in GTEX expression table

GTEx_sample_tissue_subtype : list corresponding to tissue subtypes for each sample in GTEX expression table


#### Gene-level analysis
benign.TPM : matrix (rows genes, columns subject) of benign prostate samples gene TPM

PCa.TPM  : matrix (rows genes, columns subject) of PCa samples gene TPM

CRPC.TPM  : matrix (rows genes, columns subject) of CRPC samples gene TPM

NEPC.TPM  : matrix (rows genes, columns subject) of NEPC samples gene TPM

GTEx.TPM  : matrix (rows genes, columns subject) of GTEx tissue samples gene TPM

GTEx.benign.TPM  : matrix (rows genes, columns subject) of GTEx normal prostate samples gene TPM


#### Transcript-level analysis
benign.FPKM : matrix (rows genes, columns subject) of benign prostate samples isoform FPKM

PCa.FPKM : matrix (rows genes, columns subject) of PCa samples isoform FPKM 

CRPC.FPKM : matrix (rows genes, columns subject) of CRPC samples isoform FPKM 

NEPC.FPKM : matrix (rows genes, columns subject) of NEPC samples isoform FPKM 

H660.FPKM :  matrix (rows genes, columns subject) of H660 cell line isoform FPKM

GTEX.isoform : matrix (rows genes, columns subject) of GTEx tissue samples isoform FPKM 
