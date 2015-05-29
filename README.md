# Replication material for: 'Regularized Estimation of Structural Instability in Factor Models: The US Macroeconomy and the Great Moderation'.
## Laurent Callot and Johannes Tang Kristensen.


[Link to the paper](http://lcallot.github.io/papers/ptv-fac/)

---

Date 01/06/2015

This repository contains the material necessary to replicate the empirical __application__ in: 'Regularized Estimation of Structural Instability in Factor Models: The US Macroeconomy and the Great Moderation'. 

### Required packages 

The __parsimonious__ package from the eponymous repository is required for the estimation and the __macrods__ package
provides the data:

```r
library('devtools')
install_github('lcallot/parsimonious')
library('parsimonious')
install_github('johannestang/macrods')
library('macrods')
```


In addition the following packages should be installed : _ggplot2_, _reshape2_, _xtable_, , and _zoo_, all available from CRAN.   


### /application

+ __FactorApp.R__ estimates the models and produces all output. 
+ __Libs.R__ contains helper functions. 

