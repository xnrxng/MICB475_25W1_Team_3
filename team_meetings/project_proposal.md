# Project Proposal Ideas

Idea 1: Houria!
- Dataset: Multiple Sclerosis (MS) Dataset (student-sourced!!)
- Possible research questions: (Just some very rough ideas)
- 1) Do disease-modifying therapies for MS influence the gut microbiome composition? (Compare microbiome profiles of MS patients on different medications vs. untreated individuals)
  2) How do lifestyle or dietary factors modulate the microbiome in MS patients vs. healthy controls? (Use dietary metadata or lifestyle data (e.g., fiber intake, probiotic use, smoking) and assess interactions)

Idea 2: Houria!
- Dataset: Depression Dataset (student-sourced!!)
- Possible research questions:
- 1) How does HIV/HCV coinfection alone (independent of depression) impact gut microbiome diversity and composition? (Understanding the baseline impact of coinfection on the gut microbiome is important before layering mental health on top. Compare microbiomes of HIV-only vs. HIV/HCV coinfected individuals (if the dataset allows) or vs. healthy controls)
  2) Are there identifiable microbiome-based subtypes of depression in HIV/HCV coinfected individuals? (Depression is heterogeneous, maybe some forms are more microbiome-linked than others. Cluster participants based on microbiome profiles and compare depressive symptoms or severity)
 
Idea 3: Brooke
- Dataset: Multiple Sclerosis Dataset
- Possible research questions:
- 1) How does MS treatment affect the gut microbiome? (Treatments tend to control/reduce immune system inflammation) (Basically the same as Houria's question)
  2) Do untreated patients with RRMS have different microbiomes compared to those with PPMS? (RRMS is associated with more inflammation than PPMS). How does treatment affect the microbiome in both of these cases? (I would assume that RRMS would show a greater difference because the abnormal inflammation is reduced)
 
Idea 4: Brooke
- Dataset: N/A (This is just a personal interest, we would need to find datasets if anyone is interested)
- Background:
-- Vitamin K supplementation may play a role in preventing/slowing down atherosclerosis
-- It has been suggested the Vitamin K plays a role in activating proteins that inhibit calcium deposits in tissues
-- There may be a synergistic effect with Vitamin K and Vitamin D
-- Some studies have linked Vitamin K1 intake to a lower risk of aortic stenosis, but Vitamin K2 study did not find that it was able to slow the progression
- Possible research questions:
- 1) How does vitamin K supplementation affect the gut microbiome in patients with atherosclerosis? Does it become more similar to a healthy control?
  2) Do either Vitamin K1 or Vitamin K2 have an impact on the gut microbiome and if so, are they different from each other?
  3) How does Vitamin K affect the microbiome in patients with atherosclerosis, does this change when Vitamin D is supplemented as well? How does Vitamin D alone affect the microbiome in these patients?
 
Ideas for Methods (Rui):
- For context, most of these datasets are about comparing the bacterial species of samples with different metadata/conditions. Since the format is pretty much the same, most methods can be applied on any dataset. I have done everything except networks and predicting functionality, so those two are the only ones I'm not 100% sure about.
- Basic EDA (heatmaps, scatterplots, bar charts...)
- Alpha and beta diversity analyses (boxplots, PCoa, UMAP): is there a difference in diversity between conditions?
- Differential abundance testing (I think DESeq2 can be used): which bacteria differ between conditions statistically?
- Clustering (k-means possibly): are there any subgroups of bacteria within conditions?
- Networks: are there any bacteria that tend to appear together? Is it possible to determine how they are connected?
- Predict functionality (PICRUSt): what is the functionality of this family of bacteria? Is there a pathway associated?
- Machine learning (MLR, classifiers): how well can we predict a certain attribute? Which attributes influenced the target the most?
