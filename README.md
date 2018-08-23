# Sample overlap in the context of GWAS

How much does sample overlap induce bias for different instrument strengths?

## Sim1

**Vary:**

- Confounding values
- Instrument strength
- Sample overlap

**Compare:**

- Observational estimates
- MR estimates (Wald ratios)

**Scenarios 1:**

- GWAS discovery: sample1
- Exposure: sample1
- Outcome: sample1

**Scenario 2:**

- GWAS discovery: sample1
- Exposure: sample1
- Outcome: sample2

**Scenario 3:**

- GWAS discovery: sample1
- Exposure: sample2
- Outcome: sample2

**Scenario 4:**

- GWAS discovery: sample1
- Exposure: sample2
- Outcome: sample3


Scenarios 1 and 2 are extreme versions of 3 and 4, in which winner's curse will be exaggerated.


## Sim2

Same as sim1 but more refined set of parameters

## Sim3 

Sim2 is giving similar results to those in Burgess 2016 but much noisier. Try just repeating the same thing as in Burgess 2016 - follow simulations here [https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21998](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21998) except for the following:

1. Larger range of gx effects to simulate F values more in line with GWAS results
2. Compare the use of all instruments with just those that have discovery p < 5e-8

