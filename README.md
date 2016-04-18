# PolyA

## 1. Extract annotation from organism
* Identify polyA intervals (pAi) from a given genome
* Construct array holding locations and lengths of all possible polyA events (pAi + tails)

## 2. Implement the model
* Compute distances to all pAis from BAM file
* Assign probabilities on pAis
* Estimate tail length

## 3. Tests
* Enrich unit tests

## 4. Extensions
* Introduce error bars, confidence intervals
* Implement discovery of polyA sites

---
https://www.python.org/dev/peps/pep-0008/
