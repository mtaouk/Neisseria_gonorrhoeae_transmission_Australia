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
[stats_data.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/Odds_ratio/stats_data.csv)
as input.

The input spreadsheet lists persistence in a binary system where 1 is
persistent and 0 is non persistent.

For the phenotypic antimicrobial susceptibility profiles, each
breakpoint has been assigned a value from 1 to 11 based on the
antibiotic.

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
   term          estimate std.error statistic  p.value conf.low conf.high
   <chr>            <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
 1 (Intercept)    138.      1.62     9.28     0.00232    5.80    3293.   
 2 SexM             0.279   0.347   13.5      0.000240   0.141      0.552
 3 SexOther         0.117   0.722    8.85     0.00294    0.0284     0.481
 4 AgeGroup20-29    0.912   0.345    0.0708   0.790      0.464      1.79 
 5 AgeGroup30-39    0.772   0.323    0.639    0.424      0.410      1.46 
 6 AgeGroup40-49    0.636   0.419    1.17     0.280      0.279      1.45 
 7 AgeGroup50-59    0.626   0.500    0.879    0.348      0.235      1.67 
 8 AgeGroup60-69    0.578   0.610    0.807    0.369      0.175      1.91 
 9 AgeGroup70-79    1.02    0.699    0.000727 0.978      0.259      4.01 
10 Size             1.02    0.00691  6.77     0.00925    1.00       1.03 
11 PEN              0.782   0.376    0.427    0.514      0.374      1.63 
12 TET              0.809   0.298    0.504    0.478      0.451      1.45 
13 CTRIX            0.945   0.630    0.00814  0.928      0.275      3.25 
14 CIPRO            0.866   0.162    0.785    0.376      0.630      1.19 
15 AZITH            0.383   0.354    7.37     0.00664    0.191      0.766
```

We see that sex and size of cluster are still associated with
persistence of clusters.

# GEE odds ratio model where intermediate isolates are grouped with resistant isolates

The following specified variables were included in the models: age
group, sex, size of transmission cluster, phenotypic resistance to
penicillin, phenotypic resistance to tetracycline, phenotypic resistance
to ciprofloxacin, phenotypic resistance or decreased susceptibility to
ceftriaxone and phenotypic resistance to azithromycin. For sex, an
‘unknown’ category was included to accommodate missing data, with all
other categories having complete data. In this supplementary analysis
isolates were grouped binarily as either phenotypically resistant/less
susceptible/decreased susceptibility or susceptible.

The following R code was used with
[stats_data2.csv](https://github.com/mtaouk/Neisseria_gonorrhoeae_transmission_Australia/blob/main/Supplementary_analyses/Odds_ratio/stats_data2.csv)
as input.

The input spreadsheet lists persistence in a binary system where 1 is
persistent and 0 is non persistent.

## R code:

```         
Metadata2 = read.table("stats_data2.csv", header = TRUE, sep= ",")

fit3 <- geeglm(persistant ~ Sex+AgeGroup+Size+PEN+TET+CTRIX+CIPRO+AZITH, data = Metadata2, id = Cluster, family = binomial, corstr = "independence") 

summary(fit3) 

broom::tidy(x = fit3, exp=T, conf.int = T)
```

## Results:

```         
 Coefficients:
                 Estimate Std.err  Wald Pr(>|W|)    
(Intercept)       -3.8487  1.8957  4.12   0.0423 *  
SexM              -1.0740  0.5143  4.36   0.0368 *  
SexOther/Unknown  -1.8078  0.7005  6.66   0.0099 ** 
AgeGroup20-29     -0.3327  0.3757  0.78   0.3758    
AgeGroup30-39     -0.3826  0.3364  1.29   0.2554    
AgeGroup40-49     -0.5965  0.4673  1.63   0.2018    
AgeGroup50-59     -0.7065  0.4866  2.11   0.1465    
AgeGroup60-69     -0.0922  0.6831  0.02   0.8926    
AgeGroup70-79     -0.3989  1.3806  0.08   0.7726    
Size               0.0136  0.0060  5.15   0.0232 *  
PENSUS             2.5868  0.9328  7.69   0.0056 ** 
TETSUS             0.8578  0.6361  1.82   0.1775    
CTRIXSUS          -1.7719  1.1047  2.57   0.1087    
CIPROSUS           1.6711  1.0114  2.73   0.0985 .  
AZITHSUS           4.6184  0.9599 23.15  1.5e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation structure = independence 
Estimated Scale Parameters:

            Estimate Std.err
(Intercept)    0.622   0.513
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

Here we see that sex, size of cluster and azithromycin susceptibility
are still associated with persistence of clusters. Additionally,
susceptibility to penicllin is now associated with persistent clusters.
The trends in phenotypic AMR patterns are consistent across the two
models, however we chose to use the binary model as it reflects clinical
breakpoints and is simpler to interpret.
