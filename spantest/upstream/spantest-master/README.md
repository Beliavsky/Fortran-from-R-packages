# spantest: Mean-Variance Spanning Tests in R

`spantest` is an R package providing a comprehensive suite of portfolio spanning tests for asset pricing. 
In particular, it implements classical tests by [Huberman and Kandel (1987)](https://doi.org/10.1111/j.1540-6261.1987.tb03917.x), 
[Gibbons, Ross, and Shanken (1989)](https://doi.org/10.2307/1913625), [Kempf and Memmel (2006)](https://doi.org/10.1007/BF03396737), 
[Pesaran and Yamagata (2024)](https://doi.org/10.1093/jjfinec/nbad002), and [Gungor and Luger (2016)](https://doi.org/10.1080/07350015.2015.1019510), 
alongside the recent residual-based subseries spanning tests developed by [Ardia and Sessinou (2025)](https://arxiv.org/pdf/2403.17127).

`spantest` is useful for researchers and practitioners aiming to test portfolio spanning hypotheses accounting for complex dependence structures.

Please cite the package in publications!\
By using `spantest` you agree to the following rules:

-   You must cite the package as Ardia and Seguin (2025). "spantest: Mean-Variance Spanning Tests in R", R package.
-   You assume all risk for the use of `spantest`.

Other references\
Ardia, D., & Sessinou, M. (2025).  
[Robust Inference in Large Panels and Markowitz Portfolios](https://dx.doi.org/10.2139/ssrn.5033399).  
*Working paper*.

Britten-Jones, M. (1999).\
[The Sampling Error in Estimates of Mean-Variance Efficient Portfolio Weights](https://www.jstor.org/stable/2697722).\
*Journal of Finance*, 54(2), 655-671.

Gibbons, M. R., Ross, S. A., & Shanken, J. (1989).\
[A test of the efficiency of a given portfolio](https://doi.org/10.2307/1913625).\
*Econometrica*, 57(5), 1121-1152.

Gungor, S., & Luger, R. (2016).\
[Multivariate Tests of Mean-Variance Efficiency and Spanning With a Large Number of Assets and Time-Varying Covariances](https://doi.org/10.1080/07350015.2015.1019510).\
*Journal of Business & Economic Statistics*, 34(2), 161-175.

Huberman, G., & Kandel, S. (1987).\
[Mean-variance spanning](https://doi.org/10.1111/j.1540-6261.1987.tb03917.x).\
*Journal of Finance*, 42(4), 873-888.

Kan, R., & Zhou, G. (2012).\
[Tests of Mean-Variance Spanning](https://www-2.rotman.utoronto.ca/~kan/papers/span_AEF.pdf).\
*Annals of Economics and Finance*, 13(1), 145-193.

Kempf, A., & Memmel, C. (2006).\
[Estimating the Global Minimum Variance Portfolio](https://doi.org/10.1007/BF03396737).\
*Schmalenbach Business Review*, 58, 332-348.

Pesaran, M. H., & Yamagata, T. (2024).\
[Testing for alpha in linear factor pricing models with a large number of securities.](https://doi.org/10.1093/jjfinec/nbad002).\
*Journal of Financial Econometrics*, 22(2), 407-460.

------------------------------------------------------------------------

You can install the development version of `spantest` from GitHub:

``` r
# install.packages("devtools") # if needed
devtools::install_github("ArdiaD/spantest")
```
