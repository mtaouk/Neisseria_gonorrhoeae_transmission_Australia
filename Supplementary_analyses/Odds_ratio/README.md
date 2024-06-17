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

fit2 <- geeglm(persistant ~ Sex+AgeGroup+Size+PEN+TET+CTRIX+CIPRO+AZITH, data = Metadata, id = Cluster, family = binomial, corstr = "independence") 

summary(fit2) 

broom::tidy(x = fit2, exp=T, conf.int = T)

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

```
# A tibble: 15 × 7
   term             estimate std.error statistic    p.value  conf.low conf.high
   <chr>               <dbl>     <dbl>     <dbl>      <dbl>     <dbl>     <dbl>
 1 (Intercept)        0.0213   1.90       4.12   0.0423      0.000519     0.875
 2 SexM               0.342    0.514      4.36   0.0368      0.125        0.936
 3 SexOther/Unknown   0.164    0.700      6.66   0.00986     0.0416       0.647
 4 AgeGroup20-29      0.717    0.376      0.784  0.376       0.343        1.50 
 5 AgeGroup30-39      0.682    0.336      1.29   0.255       0.353        1.32 
 6 AgeGroup40-49      0.551    0.467      1.63   0.202       0.220        1.38 
 7 AgeGroup50-59      0.493    0.487      2.11   0.147       0.190        1.28 
 8 AgeGroup60-69      0.912    0.683      0.0182 0.893       0.239        3.48 
 9 AgeGroup70-79      0.671    1.38       0.0835 0.773       0.0448      10.0  
10 Size               1.01     0.00600    5.15   0.0232      1.00         1.03 
11 PENSUS            13.3      0.933      7.69   0.00555     2.13        82.7  
12 TETSUS             2.36     0.636      1.82   0.178       0.678        8.20 
13 CTRIXSUS           0.170    1.10       2.57   0.109       0.0195       1.48 
14 CIPROSUS           5.32     1.01       2.73   0.0985      0.733       38.6  
15 AZITHSUS         101.       0.960     23.1    0.00000150 15.4        665.

```

We see that sex and size of cluster are still associated with
persistence of clusters. The trends in phenotypic AMR patterns are
consistent across the two models, however we chose to use the binary
model as it reflects clinical breakpoints and is simpler to interpret.
