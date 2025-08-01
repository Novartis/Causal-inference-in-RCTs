# Causal inference in RCTs

This repository contains code examples for several methods in a 
Causal Inference in RCTs short course. 
Novartis associates and external collaborators presented the
short course at the following conferences:

- [ICSA 2023 Applied Statistics Symposium](https://www.icsa.org/icsa-2023-applied-statistics-symposium/)
- [Joint Statistical Meetings 2023](https://ww2.amstat.org/meetings/jsm/2023/)
- [ASA Biopharmaceutical Section Regulatory-Industry Statistics Workshop 2023](https://ww2.amstat.org/meetings/biop/2023/)
- [CEN 2023](https://cen2023.github.io/home/)
- [Australian Pharmaceutical Biostatistics Group 2024](https://apbg.org.au/apbg-events/)
- [International Society for Biopharmaceutical Statistics 2024](https://www.isbiostat.org/)
- [PSI 2024](https://www.psiweb.org/conferences/about-the-conference)
- [Joint Statistical Meetings 2024](https://ww2.amstat.org/meetings/jsm/2024/)
- [International Biometric Conference 2024](https://www.ibc2024.org/home)
- [Joint Statistical Meetings 2025](https://ww2.amstat.org/meetings/jsm/2025/)

For data privacy reasons, 
the numerical results in the [`hypothetical_estimand`](hypothetical_estimand) folder are based on simulated toy datasets and will not
match the results from the short courses. 
The numerical results in the [`heart_transplant`](heart_transplant) and [`conditional_marginal`](conditional_marginal) folders
will match the results from the short courses.

# Repository contents

## Conditional and marginal effects ([`conditional_marginal`](conditional_marginal))

This folder contains example code for the 
"conditional and marginal treatment effect" lecture of the short course.

In this repository, we implemented the following approaches:

  * Conditional treatment effect point estimates and SEs using Huber-White 
  robust "sandwich" estimator. 
  * Marginal treatment effect point estimates and SEs using 
  [Ye et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10665030/)
  semiparametric approaches implemented in 
  [RobinCar2 package](https://cran.r-project.org/package=RobinCar2)
  * Functions from previous lectures, available in 
  [conditional_marginal/funs/old_funs.R](conditional_marginal/funs/old_funs.R) estimate the 
  SE of the marginal treatment effect via the following approaches:
    
      * Nonparametric bootstrap method ([Efron and Tibshirani, 1994](https://www.taylorfrancis.com/books/mono/10.1201/9780429246593/introduction-bootstrap-bradley-efron-tibshirani)) 
      * Delta method
      * Parametric bootstrap method ([Aalen et al., 1998](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0258(19971015)16:19%3C2191::AID-SIM645%3E3.0.CO;2-5)).

### How to run the scripts

To run the demo, first run [conditional_marginal/src/01_gen_data.R](conditional_marginal/src/01_gen_data.R), 
which will generate the toy dataset using the `benchtm` package.
The toy dataset has 500 samples randomized to placebo (`0`) or treatment (`1`) arm with 10 covariates. The binary response is generated from the model 
`logit(p) = 1*(X1=='Y') + 0.3*X2 + 0.3*trt`. Therefore, there are two prognostic covariates, `X1` and `X2`. This script saves the generated data to
[conditional_marginal/data/toy_data.rds](conditional_marginal/data/toy_data.rds). 
The data are also stored in the repository to ensure reproducibility. 

Second, [conditional_marginal/src/02_analysis.R](conditional_marginal/src/02_analysis.R) estimates the 
marginal and conditional treatment effects on both the risk difference and 
odds ratio scales. This script compares treatment effects under the following adjustment models:

* unadjusted, `Y ~ trt`

* adjusted with one prognostic factor `Y ~ trt + X1`

* adjusted with two prognostic factors `Y ~ trt + X1 + X2`.

## Heart transplant example ([`heart_transplant`](heart_transplant))

This folder contains two approaches to estimate 
the average causal effect on the risk difference scale (E[Y(1) - Y(0)]) for
a binary treatment Z, a binary covariate X, and a binary outcome Y.
The data example in these scripts is modified from the heart transplant example
in Chapter 1 of 
[Hernán and Robins (2020)](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/).

- [`gcomp.R`](heart_transplant/gcomp.R): Example of G-computation to estimate
E[Y(1) - Y(0)] and bootstrapping to construct a confidence interval.

- [`ipw.R`](heart_transplant/ipw.R): Example of inverse probability weighting
to estimate E[Y(1) - Y(0)] and bootstrapping to construct a confidence 
interval.

## Hypothetical estimand example ([`hypothetical_estimand`](hypothetical_estimand))

This folder contains two approaches to
estimate a hypothetical estimand. Suppose Y is an outcome, Z_0 indicates 
initial treatment assignment, and Z_1 indicates a switch to rescue medication.
Let Y(z_0, z_1) represent the potential outcome under treatment assignment
Z_0 = z_0 and rescue medication use indicated by Z_1 = z_1.
These examples estimate E[Y(1,0) - Y(0,0)], which represents the average
treatment effect in a hypothetical trial without the possibility of switching 
to rescue medication. This approach uses methods from
[Parra, Daniel, and Bartlett (2022)](https://www.tandfonline.com/doi/full/10.1080/19466315.2022.2081599).

- [`hypothetical_gcomp.R`](hypothetical_estimand/hypothetical_gcomp.R): 
Example of G-computation to estimate E[Y(1,0) - Y(0,0)] and bootstrapping to 
construct a confidence interval.

- [`hypothetical_ipw.R`](hypothetical_estimand/hypothetical_ipw.R):
Example of inverse probability weighting to estimate E[Y(1,0) - Y(0,0)] and 
bootstrapping to construct a confidence interval.

# Required packages

We list all packages that are required to run scripts within this repository.
Unless otherwise specified, packages can be installed from CRAN by using
`install.packages()`.

- [`benchtm`](https://github.com/Sophie-Sun/benchtm): 
Install by running `devtools::install_github("Sophie-Sun/benchtm")`.
- [`data.table`](https://cran.r-project.org/web/packages/data.table/index.html)
- [`future.apply`](https://cran.r-project.org/web/packages/future.apply/index.html)
- [`lmtest`](https://cran.r-project.org/web/packages/lmtest/index.html)
- [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html)
- [`progress`](https://cran.r-project.org/web/packages/progress/index.html)
- [`sandwich`](https://cran.r-project.org/web/packages/sandwich/index.html)
- [`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html)

# External links

- [Markdown document](https://oncoestimand.github.io/princ_strat_drug_dev/princ_strat_example.html) with code examples for the principal stratum estimation approaches from [Bornkamp et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1002/pst.2104). Created by Björn Bornkamp and Kaspar Rufibach.

# Code authors

- Robin Dunn, robin.dunn@novartis.com
- Jiarui Lu, jiarui.lu@novartis.com
- Tianmeng Lyu, tianmeng.lyu@novartis.com
- Tobias Muetze, tobias.muetze@novartis.com
- Cong Zhang, cong.zhang@novartis.com
