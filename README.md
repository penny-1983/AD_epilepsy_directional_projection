# AD_epilepsy_directional_projection

Analysis scripts and source data for the manuscript:

**Direction-Aware Cross-Disease Transcriptomic Analysis Reveals a Region-Sensitive Neuroimmune–Synaptic Remodeling Axis in Alzheimer’s Disease**

## Overview

This repository contains analysis scripts and source tables for a direction-aware cross-disease transcriptomic projection study in Alzheimer’s disease (AD).

The study used an epilepsy-derived directional transcriptomic reference to examine whether AD brain tissues exhibit region-sensitive molecular alignment with seizure-onset-zone-associated transcriptional states.

The analysis includes:

1. Construction of seizure onset zone (SOZ)-up and SOZ-down directional signatures from GSE268714.
2. Projection of the directional signature onto anatomically resolved AD transcriptomic data from GSE5281.
3. Calculation of sample-level EpilepsyScore using ssGSEA.
4. Region-wise and MTG-focused analyses in AD tissues.
5. Signature-size sensitivity analyses using top50, top100, and top200 directional signatures.
6. Alternative scoring robustness analysis using centered singscore.
7. Continuous association analysis between EpilepsyScore and synaptic/glia-associated module scores in MTG.
8. External cohort evaluation in GSE132903 and GSE48350.
9. Generation of main and supplementary source tables.

## Public datasets

Raw data are publicly available from the Gene Expression Omnibus (GEO):

- GSE268714
- GSE5281
- GSE132903
- GSE48350

Raw GEO files are not redistributed in this repository. Users should download the raw data directly from GEO and place them in the corresponding local data folder.

## Software

Analyses were performed in R 4.3.1.

Main R packages include:

- edgeR
- limma
- GSVA
- singscore
- fgsea
- msigdbr
- ggplot2
- pheatmap
- STRINGdb
- igraph
- openxlsx

## Repository structure

```text
scripts/
  00_setup.R
  01_build_GSE268714_directional_signature.R
  02_project_signature_GSE5281.R
  03_GSE5281_region_MTG_analysis.R
  04_sensitivity_and_singscore_robustness.R
  05_external_cohort_evaluation.R
  06_MTG_continuous_synapse_glia_analysis.R
  07_generate_figures_and_supplementary_tables.R

data/
  directional_signature/

results/
  source_tables/
  figure_source_data/
