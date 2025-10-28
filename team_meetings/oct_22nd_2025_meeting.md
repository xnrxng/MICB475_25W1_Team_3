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
- Split groups into high and low fibre with 20 g as the cutoff. Needs justification (men need x and women need y at a certain age, etc).
- Do the logistic regression first (might do RSLM? instead where CV status would be converted to 0 and 1). Then, do the LASSO at the metadata level. At the end, do the LASSO but with microbial diversity. 
- Include RDA results to project proposal as a justification of why we chose fibre and exercise. Discuss how to look at fibre in the context of exercise. Main focus should be: if you don't eat fibre but you exercise, does that rescue cardiovascular status, and viceverse? How does it change the gut microbiome?
- If fibre and exercise do not end up being significant, can explore the other variables that turned out to be significant if there is time.
- Beta diversity: split cohort into adeqante and inadequate fibre intake, and then into exercise. There should be 4 groups. Do PCoA on that. Look into the literature into what MET distribution is in order to set a good cutoff for the groupings (below a thousand and above a thousand). Do alpha diversity as well.
- DEA: check if the significant variables in the RDA have an equal distribution (if it's equal/uniform, no need to add it). Include the rest of the significant variables as covariates. Model can be run with multiple comparisons. Can put everything into one model. For instance, do those four groups and see how everything compares to one of the other groups.
- Point 3 can be ignored. Point 4 move to the top and do the LASSO with OTU table last.

**To do:**
- Look at the number of samples that would fit into every group.
- Finish proposal.
- Prepare meeting agenda for October 29th.