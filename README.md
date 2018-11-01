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



## Notes from George

i Gib/Neil
 
Thanks for this, really interesting.  A few questions.  In the summary you talk about “selecting weak instruments”, but is it not the effect size of the variants that matters, rather than they being formally weak instruments – massive sample sizes but very small variance can of course be perfectly strong instruments.

You say that “best practice would be to select SNPs that have large effects in replication”, but would it not be good to simply use the effect size from the replication?  
 
You present the distribution of F statistics from instruments for body mass index – where do these come from?  And you say you take F statistics from MR-Base, do you always have the data required to calculate these from summary GWAS data (and if so what sample that you apply them to would that be valid for? ….)? 
 
If Winner’s curse applies then you over estimate effect of the variant on the exposure and you under-estimate the causal effect of the exposure on the outcome, which isn’t stated.
 
Under further thoughts you discuss that “the expectation is that heterogeneity will be smaller using replicated effect sizes”; would this suggest that the Bowden I2 is biased upwards?  i.e. makes things look better than they are. 
 
The ratio of the “true phenotype” effect to pleiotropy magnitude would I guess be on average lowest for winner’s curse SNPs.
 
What do you intend to do with this in terms of taking forward?
 
All best

George
 