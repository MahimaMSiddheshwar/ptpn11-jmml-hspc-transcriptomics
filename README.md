# JMML vs Normal HSPCs — Bulk RNA-seq Analysis

![R](https://img.shields.io/badge/R-276A3C?style=for-the-badge&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-8DA0CB?style=for-the-badge&logo=bioconductor&logoColor=white)
![RStudio](https://img.shields.io/badge/RStudio-75AADB?style=for-the-badge&logo=rstudio&logoColor=white)
![GitHub](https://img.shields.io/badge/GitHub-181717?style=for-the-badge&logo=github&logoColor=white)
![Markdown](https://img.shields.io/badge/Markdown-000000?style=for-the-badge&logo=markdown&logoColor=white)

---


## Project Title
Inflammatory Gene Expression Signatures in PTPN11-Mutated Juvenile Myelomonocytic Leukemia Hematopoietic Stem and Progenitor Cells

## Objective
Identify differentially expressed genes and dysregulated pathways between PTPN11-mutated JMML patient HSPCs and healthy donor HSPCs using publicly available RNA-seq data (GEO: GSE183252) using the R programming Language.

## Dataset
- **Source:** GEO accession GSE183252
- **Samples:** 36 sorted HSPC populations (19 JMML, 17 Normal)
- **Original study:** Solman et al., *eLife* 2022 (https://elifesciences.org/articles/73040)

## Tools and Packages
| Purpose | Package |
|---|---|
| Data download | GEOquery |
| Differential expression | DESeq2, apeglm |
| Gene annotation | org.Hs.eg.db, EnsDb.Hsapiens.v86 |
| Pathway analysis | clusterProfiler, ReactomePA, DOSE, msigdbr |
| Visualization | ggplot2, ggrepel, enrichplot |
| Data wrangling | dplyr, tidyr, tibble |

## Pipeline Overview
1. Download counts and metadata from GEO
2. Assign disease labels (JMML vs Normal)
3. Pre-filter low-count genes (≥10 reads in ≥15 samples)
4. Run DESeq2 differential expression
5. Apply apeglm log fold-change shrinkage
6. Volcano plot of DEGs
7. KEGG over-representation analysis
8. MSigDB Hallmark GSEA
9. Reactome pathway analysis (ORA + GSEA)
10. Disease Ontology GSEA

## Thresholds
- Significance: q-value (padj) < 0.05
- Effect size: |fold change| > 1.2 (|log2FC| > 0.263)

## DEG Summary

```
 Downregulated Not Significant     Upregulated 
           2553           12259            2301 
```

## Pathway Analyses — Purpose 

- **KEGG Pathway:** Tests if DEGs over-represent known biological pathways and disease modules.
- **Hallmark GSEA:** Detects coordinated up/downregulation of curated biological processes.
- **Reactome Pathway:** Identifies enriched fine-grained signaling and regulatory pathways.
- **DOSE GSEA:** Links the DEG signature to known human disease gene patterns.

## Output Files

| File | Description |
|---|---|
| `DEGs_all.csv` | All genes with DE statistics |
| `DEGs_significant.csv` | Significant DEGs only |
| `volcano.png` | Volcano plot |
| `KEGG_dotplot.png` | KEGG enriched pathways |
| `KEGG_results.csv` | KEGG ORA result table |
| `Hallmark_barplot.png` | Hallmark GSEA NES barplot |
| `Hallmark_results.csv` | Hallmark GSEA result table |
| `Reactome_dotplot.png` | Reactome enriched pathways |
| `Reactome_results.csv` | Reactome ORA results |
| `Reactome_GSEA_results.csv` | Reactome GSEA results |
| `DOSE_GSEA_barplot.png` | Disease Ontology GSEA NES barplot |
| `DOSE_GSEA_results.csv` | Disease Ontology GSEA results |

## Plot Interpretations

**Volcano plot (`volcano.png`)**
2,301 genes up and 2,553 down in JMML; top hits TNF, NFKB1, NFKBIA, CXCL8, PTX3, MARCKSL1, ROR2 confirm strong NF-κB-driven inflammatory program in JMML HSPCs.

**Hallmark GSEA barplot (`Hallmark_barplot.png`)**
All top pathways upregulated in JMML, led by TNFα/NF-κB, interferon α/γ, and inflammatory response; (reference Solman et al.)

**KEGG dotplot (`KEGG_dotplot.png`)**
Top hits include infection pathways (Salmonella, Influenza, EBV, CMV) and metabolic stress reflecting innate immune activation; KEGG groups inflammatory genes under infection categories due to shared TLR/NF-κB gene sets.

**Reactome dotplot (`Reactome_dotplot.png`)**
Signaling by Interleukins and Diseases of growth factor receptor signaling top the list, directly mapping to JMML biology; IL-1 family and interferon signaling confirm cytokine-driven inflammation seen in patients.

**DOSE GSEA barplot (`DOSE_GSEA_barplot.png`)**
All top terms upregulated and dominated by infectious/inflammatory diseases (bacterial, viral, COVID, pneumonia, TB); JMML HSPCs transcriptionally mimic an active infection response, reinforcing the inflammatory phenotype.

## Key Findings
- 4,854 significant DEGs identified between JMML and Normal HSPCs (2,301 up, 2,553 down)
- TNF, NFKB1, NFKBIA, MARCKSL1, CXCL8, PTX3, ROR2 among top upregulated genes — classic JMML inflammatory markers
- Dominant pathway signature: TNFα via NF-κB, interferon α/γ response, inflammatory response — ref Solman et al. 2022
- JMML-defining signaling axes confirmed: JAK-STAT (IL2_STAT5, IL6_JAK_STAT3) and KRAS signaling
- Disease Ontology confirms JMML HSPCs resemble infection/inflammatory disease states transcriptionally
- Inflammatory transcriptional program reproducibly detected across DESeq2 (this study) and limma+voom (original paper)

## How to Reproduce
1. Install the packages listed under "Tools and Packages"
2. Set the working directory to the project folder
3. Run the pipeline script top to bottom — all results write to `results/GSE183252_DESeq2/`

## Reference
Solman M, *et al.* Inflammatory response in hematopoietic stem and progenitor cells triggered by activating SHP2 mutations evokes blood defects. *eLife* 2022;11:e73040. https://elifesciences.org/articles/73040
