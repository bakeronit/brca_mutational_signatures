
![made-with-R](https://img.shields.io/badge/Made%20with-R-1f425f.svg)
![GitHub last
commit](https://img.shields.io/github/last-commit/bakeronit/brca_mutational_signatures)
[![Docker
Image](https://img.shields.io/badge/docker-image-blue)](https://hub.docker.com/r/bakeronit/rstudio_hpc_r4.4)
[![CC0](https://licensebuttons.net/p/zero/1.0/88x31.png)](https://creativecommons.org/publicdomain/zero/1.0/)

# Integrating breast tumor homologous recombination deficiency status to aid germline BRCA1 and BRCA2 variant classification

This repository contains the scripts, R code, and documentation used in
the study: **“Integrating breast tumor homologous recombination
deficiency (HRD) status to aid germline BRCA1 and BRCA2 variant
classification.”**

The study evaluates tumor HRD profiles from breast cancer samples to
support interpretation of germline BRCA1/2 variants. By comparing
predictive models and calculating likelihood ratios, we aim to provide
evidence-based support for variant classification in clinical genetics.

## Repository structure

### Cohort and data resources

- 📄 [Data resources](1.data_summary.md)
- 🔧 [FFPE sample signature correction](FFPE_correction.md)

### Germline variants and BRCA1/2 classfication

- 🧬 [Germline classification](6.germline_classification.md)

### HRD analysis

- ➕ [HRDsum score calculation](4.hrdsum_calculation.md)
- 🔍 [HRD prediction using HRDetect](2.hrdetect.md)
- 🎯 [HRD prediction using CHORD](3.chord.md)
- ⚖️ [HRD prediction method comparison](5.hrdprediction_comparison.md)
- 📋 [Histopathology data and HR status correlation](7.pathology.md)

To quantify how informative HRD prediction is for classifying BRCA1/2
status, we used a likelihood ratio (LR) approach. In the main analysis,
we used CHORD-predicted HRD status to calculate LRs for HRD, BRCA1 and
BRCA2 separately. The LRs were used to assign an ACMG/AMP evidence
strength category and points using Bayesian conversions.

## R Session Info

All analyses are implemented in **R (4.4.1)** with Key package
dependencies and session details are documented in the
[rsession_info](rsession_info.md) file.
