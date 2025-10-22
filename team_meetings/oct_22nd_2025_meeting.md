# October 22nd 2025 Team 3 Meeting Notes with Chad Poloni

**Agenda:**
- Review new potential research question: How do dietary fibre intake and exercise affect cardiometabolic health and does this correspond with unique microbiome profiles?
- Ask: Do we want to include adiponectin as part of the research question?
- Ask: Determining samples to keep and samples to discard.
- Review experimental aims, rationale, and methods:
  1. Compare the diversity metrics of groups with adequate and inadequate dietary fibre intake and exercise, with cardiovascular health.
      - Group fibre intake into 3 groups. Bray-Curtis PCoA of different fibre groups, faceted by cardiovascular health. PERMANOVA for statistical significance.
      - Group exercise level into 3/4 groups. Bray-Curtis PCoA of different exercise groups, faceted by cardiovascular health. PERMANOVA for statistical significance.
      - Ask: exercise level is heavily right-skewed, group into 3 or 4 groups?
      - Ask: should we do the same thing as above, but with alpha diversity instead? In the case of alpha diversity, would we want to keep exercise and fibre as numerical or convert them to categorical?
  2. Identify taxa that differentially appear among the different fibre intake groups and exercise groups.
      - Perform DESeq2 with fibre and exercise. Perform it again but with an interaction term between fibre x cardiovascular health, and exercise x cardiovascular health.
      - Compute networks (would compute 4 networks) and compare network metrics.
      - Ask: would we include covariates for possible confounders (city, for example)?
      - Ask: confirm DESeq2 models (3 models? An interaction term?)
      - Ask: do we keep exercise and fibre as numerical or convert them to categorical?
  3. Determine pathways that are enriched across different conditions.
     - Use PICRUSt2. Compare terms by the 3 models from DESEQ2.
     - Ask: confirm 3 models again.
  4. Investigate the optimal predictors for cardiovascular health to determine if fibre and exercise fall within this category.
      - Perform redundancy analysis with the OTUs at the phylum level as response variables, metadata as explanatory variables. Assess model with a permutation test.
      - Perform logistic regression with cardiovascular health as a response variable, with the metadata as explanatory variables.
      - Perform LASSO regression with cardiovascular as a response variable. Explanatory variables should be the taxa that differentially appear from aim 2 as well as the metadata.
      - Assess each model's metrics and differences in features.
      - Ask: confirm models.

**Meeting notes:**

**To do:**
