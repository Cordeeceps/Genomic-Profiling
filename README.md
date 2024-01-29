# About the Project

This project focuses on the analysis of gene expression in pediatric brain tumor tissue samples from the Children's Brain Tumor Tissue Consortium (CBTTC). The main objective is to identify differentially expressed genes between Caucasian and African American individuals within this dataset.

# Key Project Steps

**Data Retrieval and Pre-processing:**
* Data was retrieved using the UCSCXenaTools library.
* Clinical and gene expression data were matched and processed.

**Data Normalization:**
* Gene expression data was normalized using the TMM method.

**Differential Expression Analysis:**
* Limma was used to identify differentially expressed probes.
* Probes were filtered based on a false discovery rate (FDR) of 10%.

**Top Gene Analysis:**
* The top genes with the lowest adjusted p-values were selected.
* A boxplot was generated to visualize their expression across groups.
  
**Heatmap Generation:**
* A heatmap was created to display the expression patterns of the top 30 genes across samples.

**Functional Enrichment Analysis:**
* Gene Ontology (GO) terms and KEGG pathways associated with differentially expressed genes were identified using DAVID.

# Project Summary

**Dataset:** Children's Brain Tumor Tissue Consortium (CBTTC)

**Sample Information:**
* Total Number of Probes (after filtering): 58,347
* Total Number of Samples (after filtering): 796
* Total Number of Caucasian Samples: 682
* Total Number of African American Samples: 75

**Analysis Results:**
* Number of Differentially Expressed Probes (FDR ≤ 10%): 2,172
* Top 3 Genes: AL008721.2, TBC1D3L, PWP2
* **Top GO Terms or Pathways:**
* Cholesterol biosynthetic process
* Positive regulation of insulin secretion involved in cellular response to glucose stimulus
* Negative regulation of endothelial cell proliferation
* Cell redox homeostasis
* Protein-DNA complex
* Extracellular exosome
* Protein binding
  
This project provides insights into the differential gene expression between different racial groups in pediatric brain tumor samples and identifies potential genes and pathways associated with this phenotype.

