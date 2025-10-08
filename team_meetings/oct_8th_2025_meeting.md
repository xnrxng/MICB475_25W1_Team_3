# October 8th 2025 Team 3 Meeting Notes with Chad Poloni

**Agenda:**
- Discuss the main idea of the project, which is the merging of the westernization dataset and the new depression dataset.
- Discuss how to merge these datasets (normalization, meta-analysis) and the disadvantages and the benefits.
- Where is the new depression dataset and its metadata? Does the depression dataset have annotated exercise levels (specifically METS per week)?
- Research questions:
  1. How does exercise impact the microbiome and in turn, how does it relate to depression?
     - Diversity between conditions: alpha and beta diversity
     - Bacteria that differ among conditions: differential abundance testing
     - Bacteria that tend to appear together: networks
     - The difference in functionality across conditions: PICRUSt
  2. Which attributes can explain bacterial diversity? Are we able to build an ML algorithm to predict 1) bacterial diversity or 2) condition?
 
**Meeting notes:**
- Final proposal is about individual lifestyles and health, and how it affects the microbiome
- Can look at body fat percentage and see if exercise rescues the cardiometabolic status
- Can make METS categorical (light, moderate, intense)
- Previous research has shown that exercise affects SCFA-producing bacteria
- Can look at whether exercising helps with cardiovascular health and how this is connected to body fat percentage, controlling for their diets
- Can see how cardiometabolic status changes with METS and diet. Is diet or exercise better at rescuing cardiometabolic status?
- Looking into cardiovascular condition vs doing exercise or not, balanced diet vs unbalanced
- Would need to control for different attributes (smokers, age, sex, city, etc.)
- Might need to do each analysis for each city

**To do:**
- Research into what is considered a healthy diet: percentage of carbs, fibre, protein, fat, which fats are being consumed (cholesterol, TAGs), LDL, HDL, BMI, weight, calories
- Double-check what each column actually means in the metadata
- Decide on a cut-off and bin subjects by "balanced" or "unbalanced" diet
- Do the data processing with QIIME2 without filtering, consult Chad about trimming parameters
- Analyze exercise as numerical, as well as categorical, to see if there are any trends
- Analysis would have to be per city
- Possibly, a screening assay to see which attributes are significantly changing the microbiome
- Decide the meeting agenda for Oct 15th 2025
