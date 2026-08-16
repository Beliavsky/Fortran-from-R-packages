---
title: 'gkwdist: An R Package for the Generalized Kumaraswamy Distribution Family'
tags:
  - R
  - Kumaraswamy distribution
  - bounded data
  - unit interval
  - maximum likelihood
  - RcppArmadillo
authors:
  - name: José Evandeilton Lopes
    orcid: 0009-0007-5887-4084
    affiliation: "1"
affiliations:
 - name: Federal University of Paraná (UFPR), Brazil
   index: 1
citation_author: Lopes, JE
date: 07 January 2026
year: 2026
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

Bounded continuous data on the unit interval $(0,1)$ arise frequently in scientific applications. Proportions, rates, indices, and probabilities for instance naturally fall within this domain. While the Beta distribution [@ferrari2004] and the two-parameter Kumaraswamy distribution [@kumaraswamy1980; @jones2009] serve as standard choices for such data, they may lack flexibility for datasets exhibiting heavy tails, bimodality, or complex boundary behaviors. The `gkwdist` package implements the five-parameter Generalized Kumaraswamy (GKw) distribution [@carrasco2010] and its complete hierarchy of seven nested subfamilies, providing a unified, high-performance framework for modeling bounded data with varying degrees of complexity. All computational routines are implemented in **C++** via `RcppArmadillo` [@Eddelbuettel2014], ensuring numerical stability and computational efficiency for both interactive analysis and large-scale applications.

# Statement of Need

Practitioners modeling bounded responses require distributions that balance flexibility with parsimony. The GKw family addresses this through a principled hierarchical structure: analysts can begin with simple two-parameter models (Beta, Kumaraswamy) and systematically evaluate more complex specifications (three to five parameters) only when justified by the data. This approach (starting simple and adding complexity as needed) reduces overfitting risk while accommodating challenging data features such as bimodality, asymmetric tails, and boundary concentration.

The target audience includes statisticians and data scientists working with proportions, rates, or indices in fields such as hydrology (the original domain of Kumaraswamy's distribution), econometrics, ecology, health sciences, and finance. By providing exact analytical score vectors and Hessian matrices for all seven subfamilies, `gkwdist` enables efficient maximum likelihood estimation and facilitates principled model selection via likelihood ratio tests or information criteria such as AIC and BIC.

# State of the Field

Existing R packages provide partial coverage of bounded distributions. The `extraDistr` package [@extraDistr] implements basic Kumaraswamy density, distribution, and quantile functions. The `VGAM` package [@yee2015] offers Kumaraswamy regression via the `kumar()` family function. The `betareg` package [@betareg] provides comprehensive Beta regression modeling with diagnostic tools. The `gamlss.dist` package [@rigby2019] includes Beta variants within the GAMLSS framework, supporting distributional regression for location, scale, and shape parameters. However, none of these packages provides the complete GKw hierarchy covering the seven nested distributions ranging from two to five parameters with analytical derivatives essential for efficient inference. The `gkwdist` package fills this gap by offering: (i) all subfamilies in a unified API following standard R conventions; (ii) exact analytical gradients and Hessians for each distribution; and (iii) numerically stable evaluation across challenging parameter regions via robust logarithmic transformations.

## Mathematical Framework

Let $X \sim \text{GKw}(\alpha, \beta, \gamma, \delta, \lambda)$ with all parameters positive. The probability density function is:

$$f(x; \boldsymbol{\theta}) = \frac{\lambda \alpha \beta}{B(\gamma, \delta+1)} x^{\alpha-1} (1-x^\alpha)^{\beta-1} w^{\gamma\lambda-1} (1-w^\lambda)^\delta$$

where $w = 1 - (1-x^\alpha)^\beta$, $\boldsymbol{\theta} = (\alpha, \beta, \gamma, \delta, \lambda)^\top$, and $B(\cdot, \cdot)$ denotes the Beta function. The cumulative distribution function is $F(x) = I_{w^\lambda}(\gamma, \delta+1)$, where $I_y(a,b)$ represents the regularized incomplete beta function.

Seven subfamilies arise through parameter constraints (\autoref{tab:subfamilies}). The Kumaraswamy (Kw) and Exponentiated Kumaraswamy (EKw) subfamilies admit closed-form quantile functions, enabling efficient random variate generation via inversion. The Beta distribution emerges as a special case when $\alpha = \beta = \lambda = 1$, establishing a direct connection to the foundational work of @ferrari2004 with a simple variation in its second parameter $\alpha = \gamma$ and $\beta = \delta + 1$.

: GKw subfamily hierarchy with parameter constraints and quantile tractability. \label{tab:subfamilies}

| Subfamily | Parameters | Constraint | Closed Quantile |
|:----------|:----------:|:-----------|:---------------:|
| GKw (`gkw`) | 5 | — | No |
| BKw (`bkw`) | 4 | $\lambda=1$ | No |
| KKw (`kkw`) | 4 | $\gamma=1$ | Yes |
| EKw (`ekw`) | 3 | $\gamma=1, \delta=0$ | Yes |
| Mc (`mc`) | 3 | $\alpha=\beta=1$ | No |
| Kw (`kw`) | 2 | $\gamma=\delta=0, \lambda=1$ | Yes |
| Beta (`beta_`)| 2 | $\alpha=\beta=\lambda=1$ | No |

# Software Design

The package prioritizes computational efficiency and numerical stability through complete C++ implementation via `RcppArmadillo` [@Eddelbuettel2014]. Key design decisions include:

1. **Unified API**: All distribution functions (`dgkw`, `pgkw`, `qgkw`, `rgkw` and analogous functions for each subfamily) and inference functions (negative log-likelihood `llgkw`, gradient `grgkw`, Hessian `hsgkw`) follow standard R naming conventions.

2. **Analytical derivatives**: Exact score vectors and Hessian matrices enable gradient-based optimization via `optim()` with method `"BFGS"`, `"L-BFGS-B"` among others, substantially improving convergence speed compared to numerical differentiation.

3. **Numerical stability**: Logarithmic transformations handle extreme parameter values and boundary conditions, preventing underflow/overflow in density calculations.

```r
library(gkwdist)
x <- rkw(1000, alpha = 2, beta = 3)
fit <- optim(c(1, 1), fn = llkw, gr = grkw, data = x, method = "BFGS")
```

Mathematical derivations of all analytical gradients and Hessians are documented in the package vignettes available at the package website.

# Research Impact Statement

The `gkwdist` package serves as the computational foundation for `gkwreg` [@gkwreg], which implements regression models for bounded responses using the GKw family. The package has been available on CRAN since November 2025, with source code hosted on GitHub under the MIT license. The complete test suite achieves over 95% code coverage, and all analytical derivatives have been validated against numerical differentiation via the `numDeriv` package with near zero precision.

Applications span hydrology [@kumaraswamy1980], econometrics, ecology, and health sciences, domains where bounded data modeling is essential. The hierarchical structure facilitates systematic model comparison, allowing researchers to identify the most parsimonious adequate model for their specific application.

# AI Usage Disclosure

Documentation and manuscript text were refined with assistance from Claude (Anthropic, version Opus 4.5). Mathematical derivations of score vectors and Hessian matrices were independently verified by systematic comparison with numerical differentiation via the `numDeriv` package across multiple parameter configurations. The core C++ implementation represents substantial original development effort by the author.

# Acknowledgements

The author thanks Prof. Wagner Hugo Bonat (UFPR) for guidance on statistical methodology and the R community for valuable feedback during package development.

# References