# October 29th 2025 Team 3 Meeting Notes with Chad Poloni

**Agenda:**
- Project proposal draft so far: https://docs.google.com/document/d/1ChwnbNqw5gJSZl0BKtLdcWXla6qauKmJxd-17VJKh8g/edit?tab=t.0
- Groups (result artifact 8) are uneven
- RSLM instead of logistic regression?
- Include logistic regression stats? (MET p.value < 0.1)
- LASSO with OTU table at the end, part of aim 1 or a different aim?

**Meeting notes:**
- Condense the proposal. Leave aim number 4 in if we can condense the whole thing but if not, take it out. 
- Move the lasso with the OTU (last point in aim 1) to aim 3 (at the beginning, then do deseq and compare).
- Make sure you add citations.
- Logistic regression (linear model). Transform and do parametric test.
- Need to normalize the data (arc sign? or another method? Chad said he'll ask and update us?) before doing a test (like linear model) that would be done on a parametric dataset (which this obvs is not). 
- Should normalize before doing RDA and lasso and linear model.
- Package name tells you about the types of normalization.
- When doing permonova and RDA (using vegan package), it would ask you for methods (to normalize?) so make sure don’t normalize again when running that!
- We have to run a diablo if we want to merge/integrate different types of data

**To do:**
- Finish proposal
- Agenda for Nov 5th meeting
- Get started on the analysis (Houria: alpha diversity. Quin: beta diversity. Rui: predictors. Brooke: deseq)
