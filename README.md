
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
- 📋 [Pathology and HR status](7.pathology.md)

### Germline variants and BRCA1/2 classfication

- 🧬 [Germline classification](6.germline_classification.md)

### HRD analysis

- ➕ [HRDsum score calculation](4.hrdsum_calculation.md)
- 🔍 [HRD prediction using HRDetect](2.hrdetect.md)
- 🎯 [HRD prediction using CHORD](3.chord.md)
- ⚖️ [HRD prediction method comparison](5.hrdprediction_comparison.md)

To quantify the informativeness of HRD prediction for BRCA1/2
classification, we employed a likelihood ratio (LR) approach. In main
results, we used CHORD-predicted HRD status, we calculated LRs for
BRCA1/2.

## R Session Info

All analyses are implemented in **R (≥4.4.1)** with Key package
dependencies and session details are documented in the
[rsession_info](rsession_info.md) file.
