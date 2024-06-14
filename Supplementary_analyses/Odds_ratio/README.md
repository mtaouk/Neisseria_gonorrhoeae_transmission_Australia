# GEE odds ratio model using continuous MIC variables

To assess variables associated with cluster persistence, we applied a
multivariable logistic regression model via generalised estimating
equation (GEE) to calculate adjusted odds ratios using an independence
model with geepack (v1.3.9). The following specified variables were
included in the models: age group, sex, size of transmission cluster,
phenotypic resistance to penicillin, phenotypic resistance to
tetracycline, phenotypic resistance to ciprofloxacin, phenotypic
resistance or decreased susceptibility to ceftriaxone and phenotypic
resistance to azithromycin. For sex, an ‘unknown’ category was included
to accommodate missing data, with all other categories having complete
data. In this supplementary analysis we used continuous MIC values for
the phenotypic antimicrobial susceptibility profiles.

The following R code was used with
s[tats_data.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/Odds_ratio/stats_data.csv)
as input.

The input spreadsheet lists persistence in a binary system where 1 is
persistent and 0 is non persistent.

## R code:

```         
library(geepack)

Metadata = read.table("stats_data.csv", header = TRUE, sep= ",")

fit3 <- geeglm(persistant ~ Sex+AgeGroup+Size+PEN+TET+CTRIX+CIPRO+AZITH, data = Metadata, id = Cluster, family = binomial, corstr = "independence") 

summary(fit3) 
```

## Results:

```         
 Coefficients:
               Estimate   Std.err   Wald Pr(>|W|)    
(Intercept)    4.928596  1.617803  9.281  0.00232 ** 
SexM          -1.275334  0.347233 13.490  0.00024 ***
SexOther      -2.147489  0.721981  8.847  0.00294 ** 
AgeGroup20-29 -0.091852  0.345161  0.071  0.79015    
AgeGroup30-39 -0.258350  0.323177  0.639  0.42405    
AgeGroup40-49 -0.452894  0.419387  1.166  0.28019    
AgeGroup50-59 -0.468398  0.499611  0.879  0.34849    
AgeGroup60-69 -0.547614  0.609648  0.807  0.36905    
AgeGroup70-79  0.018856  0.699239  0.001  0.97849    
Size           0.017994  0.006914  6.774  0.00925 ** 
PEN           -0.245739  0.376097  0.427  0.51350    
TET           -0.211504  0.297985  0.504  0.47784    
CTRIX         -0.056818  0.629740  0.008  0.92811    
CIPRO         -0.143813  0.162365  0.785  0.37576    
AZITH         -0.959885  0.353599  7.369  0.00664 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation structure = independence 
Estimated Scale Parameters:

            Estimate Std.err
(Intercept)   0.6458    2.57
Number of clusters:   31  Maximum cluster size: 709 
```

We see that sex and size of cluster are still associated with
persistence of clusters. The trends in phenotypic AMR patterns are
consistent across the two models, however we chose to use the binary
model as it reflects clinical breakpoints and is simpler to interpret.
