# NeuroSurfacePro

Multi-step bioinformatics pipeline to identify NEPC-specific, over-expressed gene transcripts that encode surface proteins.

This work is further described in: [citation to be added upon publication]

### Data availability
RNAseq data was obtained from  a previously published (Beltran et al., 2016; Berger et al, 2019) clinical cohort of tissue from NEPC tumors (n = 27), benign prostate tissue (n = 66), locally advanced prostate cancer (PCA, n = 68), and mCRPC (n = 72).

The GTEx data used for the analyses described were obtained from the [GTEx Portal](https://gtexportal.org/home/) on 03/01/2015 and/or dbGaP accession number phs000424.vN.pN on 03/01/2015


### Gene-Level Pipeline expected input:

benign.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of benign prostate samples

PCa.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of PCa samples

CRPC.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of CRPC samples

NEPC.FPKM : matrix (rows genes, columns subject) contain isoform FPKM of NEPC samples

H660.FPKM :  matrix (rows genes, columns subject) contain isoform FPKM of H660 cell line


GTEX.isoform : matrix (rows genes, columns subject) contain isoform FPKM of GTEx tissue samples

GTEx_sample_tissue : list corresponding to tissue types for each sample in GTEX.isoform

GTEx_sample_tissue_subtype2 : list corresponding to tissue subtypes for each sample in GTEX.isoform


NEPCvBenign : isoform-level cuffdiff result data frame of NEPC compared to benign prostate samples

NEPCvPCa : isoform-level cuffdiff result data frame of NEPC compared to PCa samples

NEPCvCRPC : isoform-level cuffdiff result data frame of NEPC compared to CRPC samples


redflag_tissues <- table containing tissue importance labels for each tissue in GTEx ( 0 = dispensible, 1 = non-critical, 2 = critical )
