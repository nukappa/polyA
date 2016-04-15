# PolyA

## 1. Extract annotation from organism
* Identify polyA intervals (pAi) from a given genome
* Extract 3' UTR information from a given gtf 
* Construct array holding locations and lengths of all possible polyA events (pAi + tails)

## 2. Implement the model
* Discretize according to bioanalyzer profile
* Compute distances to all pAis from BAM file
* Assign probabilities on pAis
* Estimate tail length

## 3. Tests
* Upload test data (bioanalyzer + small set of distances)
* Write unit tests

## 4. Extensions
* Introduce error bars, confidence intervals
* Implement discovery of polyA sites