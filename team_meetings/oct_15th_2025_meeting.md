# October 15th 2025 Team 3 Meeting Notes with Chad Poloni

**Agenda:**
- Discussing our chosen "healthy diet": percentage of carbs, fibre, protein, fat, which fats are being consumed (cholesterol, TAGs), LDL, HDL, BMI, weight, (not including calories?)
- Discuss data normalization (log transform, CPM, TPM) and batch correction
- Discussing our specific research question(s). Some potential options are:
 1. How does the gut microbiota change between a balanced diet and an unbalanced diet, across different exercise levels? How is cardiovascular health affected by this?
    - Which set of species/attributes predicts cardiovascular health the best? Can cardiovascular health even be predicted? Append OTU table and/or diversity metrics to patient metadata, split data into 2, perform LASSO regression on the training dataset for factor selection. Apply to the testing dataset for model metrics. Target variable is cardiovascular health. Possibly look into other models.
    - Is there a difference in diversity? Compare the diversity metrics of both diets, as well as of different exercise groups.
    - Are there any species that differentially appear? Differential abundance testing of both diets, as well as of different exercise groups as contrasts. Also test as interaction term. See which species overlap between the two.
    - What is the functionality of this group of microbiota? Predict metagenomic functional profiles, compare functions by diet/exercise and link to CV health.
    - Are there any species that appear consistently together? Identify taxa that co-occur consistently, and whether network structure differs by diet/exercise.
   
 2. How does cardiovascular health relate to adiponectin levels and microbiome profile? How does exercise influence this?
    - Adiponectin has anti-inflammatory properties protect vascular system, heart, lungs, colon and these levels vary based on sex and BMI
    - Obesity is related to low adiponectin levels, and both obesity and adiponectin are related to poor cardiovascular health
    - Found another paper that looks into adiponectin levels, microbiome, and obesity which found reduced adiponectin in all groups with any type of metabolic abnormalities compared to healthy individuals
       - This paper did not incorporate exercise
  
- Columns of the metadata
 1. Adiponectin: a hormone from fat tissue that helps with insulin sensitivity and inflammation, and affects metabolic processes like fatty acid breakdown and regulating glucose levels.
 2. Diastolic bp: diastolic blood pressure is the bottom number in a blood pressure reading, representing the pressure in your arteries when your heart is at rest between beats. High diastolic pressure can indicate increased pressure in your arteries and may be linked to serious health complications if left unmanaged. A healthy diastolic reading is typically less than 80 mm Hg.
 3. Hemoglobin a1c: a Hemoglobin A1c (HbA1c) test is a blood test measuring your average blood sugar (glucose) over the past two to three months by checking the amount of glucose attached to hemoglobin in your red blood cells. It's used to diagnose and monitor prediabetes and diabetes, with results reported as a percentage. A normal A1c is typically below 5.7%, while 5.7% to 6.4% indicates prediabetes, and 6.5% or higher suggests diabetes.
 4. CRP: short for C-reactive protein, a protein produced by the liver that indicates the presence of inflammation in the body. A CRP test is a blood test used to measure levels of this protein to help detect and monitor infections, inflammatory diseases, and tissue damage.
 5. Systolic bp: systolic blood pressure is the top number in a blood pressure reading, representing the force of blood against artery walls when the heart contracts. Normal blood pressure is often considered to be around 120/80 mmHg. High blood pressure (hypertension) is consistently high readings, which can increase the risk of heart disease. For people over 50, a higher systolic number is a particularly important indicator of heart disease risk. 
 
**Meeting notes:**

**To do:**
