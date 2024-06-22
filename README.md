# Missing Data Corresponding Curve Regression (mdCCR)

This code is for paper titled- "Assessing Reproducibility of High-throughput Experiments in the Case of Missing Data".
Link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9039958/

## Introduction

High-throughput experiments are an essential part of modern biological and biomedical research.
The outcomes of high-throughput biological experiments often have a lot of missing observations
due to signals below detection levels. For example, most single-cell RNA-seq (scRNA-seq)
protocols experience high levels of dropout due to the small amount of starting material, leading to
a majority of reported expression levels being zero. Though missing data contain information
about reproducibility, they are often excluded in the reproducibility assessment, potentially
generating misleading assessments.

Here, we propose mdCCR, a regression model to assess how the reproducibility of high-throughput
experiments is affected by the choices of operational factors (e.g., platform or sequencing depth)
when a large number of measurements are missing. Using a latent variable approach, we extend
correspondence curve regression (CCR), a recently proposed method for assessing the effects of
operational factors to reproducibility, to incorporate missing values. Using simulations, we show
that our method is more accurate in detecting differences in reproducibility than existing measures
of reproducibility. We illustrate the usefulness of our method using a single-cell RNA-seq dataset
collected on HCT116 cells. We compare the reproducibility of different library preparation
platforms and study the effect of sequencing depth on reproducibility, thereby determining the
cost-effective sequencing depth that is required to achieve sufficient reproducibility.

## Scripts

Nelsen_mdCCR.R : This script uses mdCCR to estimate the reproducibility in case of data following Nelsen 4.2.12 copula distribution. The results correspond to section 4.3 in the paper.

Gumbel_mdCCR.R : This script uses mdCCR to estimate the reproducibility in case of data following Gumbel-Hougaard copula distribution. The results correspond to section 4.3 in the paper.

Realdata_platform_mdCCR.R : This script uses mdCCR to estimate the reproducibility of different platforms (i) TransPlex, (ii) SMARTer, and (iii) C1 to produce scRNA-seq data. The results correspond to section 5.1 in the paper.

Realdata_downsample_mdCCR.R : This script uses mdCCR to estimate the reproducibility of different sequencing depths to produce scRNA-seq data. The results correspond to section 5.2 in the paper.


