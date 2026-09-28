[![R-CMD-check](https://github.com/nilforooshan/mbend/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nilforooshan/mbend/actions/workflows/R-CMD-check.yaml)
[![cran-version-last-release](https://www.r-pkg.org/badges/version-last-release/mbend)](https://cran.r-project.org/package=mbend)
[![cran-downloads-total](https://cranlogs.r-pkg.org/badges/grand-total/mbend?color=green)](https://cran.r-project.org/package=mbend)
[![cran-downloads-month](https://cranlogs.r-pkg.org/badges/last-month/mbend?color=green)](https://cran.r-project.org/package=mbend)
[![cran-downloads-week](https://cranlogs.r-pkg.org/badges/last-week/mbend?color=green)](https://cran.r-project.org/package=mbend)
[![github-release](https://img.shields.io/github/release/nilforooshan/mbend.svg)](https://github.com/nilforooshan/mbend)
[![github-license](https://img.shields.io/github/license/nilforooshan/mbend.svg)](https://github.com/nilforooshan/mbend/blob/master/LICENSE)
[![github-contributors](https://img.shields.io/github/contributors/nilforooshan/mbend.svg)](https://github.com/nilforooshan/mbend/graphs/contributors/)

# mbend: An R package for bending non-positive-definite matrices to positive-definite

Mohammad Ali Nilforooshan

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-0339-5442" href="https://orcid.org/0000-0003-0339-5442" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0003-0339-5442</a></div>

---

## Description

_Bending non-positive-definite (symmetric) matrices to positive-definite, using weighted and unweighted methods_

The `mbend` package is used for bending symmetric non-positive-definite matrices to positive-definite (PD). For a matrix to be invertible, it has to be PD. The methods of Jorjani et al. (2003) and Schaeffer (2014) are used in this package. The unweighted method of Schaeffer (2014) is also extended to weighted bending, with the possibility of choosing between unweighted and weighted bending.

## Application

Start with loading the package library:

```r
library(mbend)
```

Consider the following non-PD covariance matrix (Jorjani et al., 2003):

```r
V = matrix(nrow=5, ncol=5, c(
  100,  95,  80,  40,  40,
   95, 100,  95,  80,  40,
   80,  95, 100,  95,  80,
   40,  80,  95, 100,  95,
   40,  40,  80,  95, 100))
```

The following command can be used to bend matrix `V` to a PD matrix:

```r
bend(V)
#> Unweighted bending
#> max.iter = 10000
#> small.positive = 1e-04
#> method = hj
#> Convergence met after 1 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 103.16918  90.83313  79.47122  44.53731  37.07254
#> [2,]  90.83313 106.49734  94.18961  74.07039  44.53731
#> [3,]  79.47122  94.18961 102.31355  94.18961  79.47122
#> [4,]  44.53731  74.07039  94.18961 106.49734  90.83313
#> [5,]  37.07254  44.53731  79.47122  90.83313 103.16918
#>
#> $init.ev
#> [1] 399.475997  98.523500  23.646897  -3.122893 -18.523500
#>
#> $final.ev
#> [1] 399.4760  98.5235  23.6469   0.0001   0.0001
#>
#> $min.dev
#> [1] -5.929612
#>
#> $max.dev
#> [1] 6.497343
#>
#> $loc.min.dev
#> row col
#>   2   4
#>
#> $loc.max.dev
#> row col
#>   4   4
#>
#> $ave.dev
#> [1] 0.7234705
#>
#> $AAD
#> [1] 3.372692
#>
#> $Cor
#> [1] 0.985632
#>
#> $RMSD
#> [1] 3.927457
```

The above command is equivalent to `bend(inmat=V)`, where `inmat` is the argument that takes the matrix to be bent, or the following command:

```r
bend(V, max.iter=10000, small.positive=0.0001, method="hj")
```

This runs the unweighted bending method of Jorjani et al. (2003) (`method="hj"`), with maximum 10000 number of iterations (`max.iter=10000`), and eigenvalues smaller than 0.0001 are replaced with this small positive value (`small.positive=0.0001`). Providing the default parameters, the corresponding arguments can be omitted (e.g., `bend(V)`).

The output object is a list of several items, listed below:

- bent : The bent `matrix`.
- init.ev : Eigenvalues of the initial (`inmat`) matrix.
- final.ev : Eigenvalues of the `bent` matrix.
- min.dev : `min(bent - inmat)`.
- max.dev : `max(bent - inmat)`.
- loc.min.dev : Location (indices) of `min.dev` element.
- loc.max.dev : Location (indices) of `max.dev` element.
- ave.dev : Average deviation (`bent - inmat`) of the upper triangle elements (excluding diagonal elements for correlation matrices).
- AAD : Average absolute deviation of the upper triangle elements (excluding diagonal elements for correlation matrices) of `bent` and `inmat`.
- Cor : Correlation between the upper triangle elements (excluding diagonal elements for correlation matrices) of `bent` and `inmat`.
- RMSD : Root of mean squared deviation of the upper triangle elements (excluding diagonal elements for correlation matrices) of `bent` and `inmat`.

There might be different precision involved with different elements of a non-PD matrix. In this case, a weighted bending is recommended. Jorjani et al. (2003) used the reciprocal of the number data points in common between pairs of variables, as weights. Considering the following matrix for the number of data points in common between variables (Jorjani et al., 2003):

```r
W = matrix(nrow=5, ncol=5, c(
  1000,  500,   20,   50,  200,
   500, 1000,  500,    5,   50,
    20,  500, 1000,   20,   20,
    50,    5,   20, 1000,  200,
   200,   50,   20,  200, 1000))
```

Matrix `V` is bent using the following command:

```r
#> bend(inmat=V, wtmat=W, reciprocal=TRUE)
#> Weighted bending
#> reciprocal = TRUE
#> max.iter = 10000
#> small.positive = 1e-04
#> method = hj
#> Convergence met after 428 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.16154  94.51747  82.93843  43.56681  39.18292
#> [2,]  94.51747 100.62488  93.98150  59.99921  45.84555
#> [3,]  82.93843  93.98150 100.69547  84.89294  73.13431
#> [4,]  43.56681  59.99921  84.89294 100.30697  94.23170
#> [5,]  39.18292  45.84555  73.13431  94.23170 100.17942
#>
#> $init.ev
#> [1] 399.475997  98.523500  23.646897  -3.122893 -18.523500
#>
#> $final.ev
#> [1] 3.876377e+02 1.020313e+02 1.229913e+01 7.667873e-05 1.210388e-06
#>
#> $min.dev
#> [1] -20.00079
#>
#> $max.dev
#> [1] 5.84555
#>
#> $loc.min.dev
#> row col
#>   2   4
#>
#> $loc.max.dev
#> row col
#>   2   5
#>
#> $ave.dev
#> [1] -1.716058
#>
#> $AAD
#> [1] 3.625268
#>
#> $Cor
#> [1] 0.9623206
#>
#> $RMSD
#> [1] 6.368693
#>
#> $w_gt_0
#> [1] 15
#>
#> $wAAD
#> [1] 0.6100103
#>
#> $wCor
#> [1] 0.9955121
#>
#> $wRMSD
#> [1] 0.5326555
```

Using `wtmat=1/W, reciprocal=FALSE`, the argument `reciprocal` could be omitted, because `FALSE` is the default parameter for the argument `reciprocal`. For the same reason, `max.iter=10000, small.positive=0.0001, method="hj"` are omitted, unless different parameters are provided to these arguments.

If there is high confidence about some elements of the non-PD matrix to remain unchanged after bending, the corresponding weights are set to zero. For example, to keep the first 2 &times; 2 block of `V` unchanged during the bending procedure:

```r
W2 = W; W2[1:2, 1:2] = 0
bend(V, W2, reciprocal=TRUE)
#> Weighted bending
#> reciprocal = TRUE
#> max.iter = 10000
#> small.positive = 1e-04
#> method = hj
#> Convergence met after 475 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.00000  95.00000  83.70895  43.70489  39.12852
#> [2,]  95.00000 100.00000  93.89940  59.58422  46.15038
#> [3,]  83.70895  93.89940 100.73015  84.47456  72.91450
#> [4,]  43.70489  59.58422  84.47456 100.31566  94.21299
#> [5,]  39.12852  46.15038  72.91450  94.21299 100.18322
#>
#> $init.ev
#> [1] 399.475997  98.523500  23.646897  -3.122893 -18.523500
#>
#> $final.ev
#> [1] 3.876582e+02 1.021544e+02 1.141629e+01 7.646818e-05 1.174802e-07
#>
#> $min.dev
#> [1] -20.41578
#>
#> $max.dev
#> [1] 6.150379
#>
#> $loc.min.dev
#> row col
#>   4   2
#>
#> $loc.max.dev
#> row col
#>   2   5
#>
#> $ave.dev
#> [1] -1.732836
#>
#> $AAD
#> [1] 3.705269
#>
#> $Cor
#> [1] 0.9597991
#>
#> $RMSD
#> [1] 6.564345
#>
#> $w_gt_0
#> [1] 12
#>
#> $wAAD
#> [1] 0.7705437
#>
#> $wCor
#> [1] 0.9951209
#>
#> $wRMSD
#> [1] 0.6080559
```

For weighted bending, we get extra statistics in the output:

- `w_gt_0` : Number of weight elements greater than 0, in the upper triangle of `wtmat` (for weighted bending).
- `wAAD` : Weighted `AAD` (for weighted bending).
- wCor : Weighted `Cor` (for weighted bending).
- wRMSD : Weighted `RMSD` (for weighted bending).

To bend `V` using the method of Schaeffer (2014):

```r
#> bend(inmat=V, method="lrs")
#> Unweighted bending
#> max.iter = 10000
#> method = lrs
#> Convergence met after 1 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 103.18978  90.82704  79.43676  44.56754  37.06769
#> [2,]  90.82704 106.54177  94.13680  74.06295  44.56754
#> [3,]  79.43676  94.13680 102.46429  94.13680  79.43676
#> [4,]  44.56754  74.06295  94.13680 106.54177  90.82704
#> [5,]  37.06769  44.56754  79.43676  90.82704 103.18978
#>
#> $init.ev
#> [1] 399.475997  98.523500  23.646897  -3.122893 -18.523500
#>
#> $final.ev
#> [1] 399.47599653  98.52349955  23.64689694   0.20358328   0.07740478
#>
#> $min.dev
#> [1] -5.937047
#>
#> $max.dev
#> [1] 6.54177
#>
#> $loc.min.dev
#> row col
#>   4   2
#>
#> $loc.max.dev
#> row col
#>   4   4
#>
#> $ave.dev
#> [1] 0.732955
#>
#> $AAD
#> [1] 3.408707
#>
#> $Cor
#> [1] 0.9854511
#>
#> $RMSD
#> [1] 3.954199
```

The method of Schaeffer (2014) does not require the argument `small.positive`, and this argument is ignored. This method is originally an unweighted bending method. However, in this package, it is extended to accommodate weighted bending (i.e., a combination of Schaeffer (2014) and Jorjani et al. (2003) methods). Weighted bending of `V` with reciprocals of `W` using the method of Schaeffer (2014):

```r
bend(V, W, reciprocal=TRUE, method="lrs")
#> Weighted bending
#> reciprocal = TRUE
#> max.iter = 10000
#> method = lrs
#> Convergence met after 787 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 100.16194  94.51554  82.97457  43.57119  39.18071
#> [2,]  94.51554 100.62754  93.97669  60.00766  45.86380
#> [3,]  82.97457  93.97669 100.69806  84.86261  73.10529
#> [4,]  43.57119  60.00766  84.86261 100.30785  94.22946
#> [5,]  39.18071  45.86380  73.10529  94.22946 100.18003
#>
#> $init.ev
#> [1] 399.475997  98.523500  23.646897  -3.122893 -18.523500
#>
#> $final.ev
#> [1] 3.876365e+02 1.020234e+02 1.228261e+01 3.242161e-02 4.435327e-04
#>
#> $min.dev
#> [1] -19.99234
#>
#> $max.dev
#> [1] 5.863797
#>
#> $loc.min.dev
#> row col
#>   4   2
#>
#> $loc.max.dev
#> row col
#>   5   2
#>
#> $ave.dev
#> [1] -1.715806
#>
#> $AAD
#> [1] 3.633803
#>
#> $Cor
#> [1] 0.9622401
#>
#> $RMSD
#> [1] 6.374766
#>
#> $w_gt_0
#> [1] 15
#>
#> $wAAD
#> [1] 0.6122026
#>
#> $wCor
#> [1] 0.9954884
#>
#> $wRMSD
#> [1] 0.534688
```

Function `bend` automatically considers any matrix with all diagonal elements equal to one, as a correlation matrix. Considering the correlation matrix `V2` (`V` converted to a correlation matrix):

```r
V2 = cov2cor(V)
bend(V2, W, reciprocal=TRUE)
#> Weighted bending
#> reciprocal = TRUE
#> max.iter = 10000
#> small.positive = 1e-04
#> method = hj
#> Found a correlation matrix.
#> Convergence met after 286 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9447536 0.8342673 0.4377103 0.3912583
#> [2,] 0.9447536 1.0000000 0.9385087 0.6004642 0.4630373
#> [3,] 0.8342673 0.9385087 1.0000000 0.8394178 0.7248761
#> [4,] 0.4377103 0.6004642 0.8394178 1.0000000 0.9418815
#> [5,] 0.3912583 0.4630373 0.7248761 0.9418815 1.0000000
#>
#> $init.ev
#> [1]  3.99475997  0.98523500  0.23646897 -0.03122893 -0.18523500
#>
#> $final.ev
#> [1] 3.868843e+00 1.015420e+00 1.156566e-01 7.740771e-05 2.031735e-06
#>
#> $min.dev
#> [1] -0.1995358
#>
#> $max.dev
#> [1] 0.06303735
#>
#> $loc.min.dev
#> row col
#>   2   4
#>
#> $loc.max.dev
#> row col
#>   5   2
#>
#> $ave.dev
#> [1] -0.02838248
#>
#> $AAD
#> [1] 0.05538548
#>
#> $Cor
#> [1] 0.9463128
#>
#> $RMSD
#> [1] 0.0803483
#>
#> $w_gt_0
#> [1] 10
#>
#> $wAAD
#> [1] 0.01416959
#>
#> $wCor
#> [1] 0.9942541
#>
#> $wRMSD
#> [1] 0.01074558
```

For correlation matrices, diagonal elements are not used in obtaining the statistics.
Because the argument `method` is not provided, the default parameter `"hj"` is used. To do the same using the method of Schaeffer (2014):

```r
bend(V2, W, reciprocal=TRUE, method="lrs")
#> Weighted bending
#> reciprocal = TRUE
#> max.iter = 10000
#> method = lrs
#> Found a correlation matrix.
#> Convergence met after 700 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9447383 0.8345316 0.4377385 0.3912447
#> [2,] 0.9447383 1.0000000 0.9384618 0.6005923 0.4631893
#> [3,] 0.8345316 0.9384618 1.0000000 0.8392102 0.7245638
#> [4,] 0.4377385 0.6005923 0.8392102 1.0000000 0.9418734
#> [5,] 0.3912447 0.4631893 0.7245638 0.9418734 1.0000000
#>
#> $init.ev
#> [1]  3.99475997  0.98523500  0.23646897 -0.03122893 -0.18523500
#>
#> $final.ev
#> [1] 3.868820e+00 1.015338e+00 1.155776e-01 2.617687e-04 2.856950e-06
#>
#> $min.dev
#> [1] -0.1994077
#>
#> $max.dev
#> [1] 0.06318928
#>
#> $loc.min.dev
#> row col
#>   4   2
#>
#> $loc.max.dev
#> row col
#>   5   2
#>
#> $ave.dev
#> [1] -0.02838562
#>
#> $AAD
#> [1] 0.0554775
#>
#> $Cor
#> [1] 0.946235
#>
#> $RMSD
#> [1] 0.08039991
#>
#> $w_gt_0
#> [1] 10
#>
#> $wAAD
#> [1] 0.01420763
#>
#> $wCor
#> [1] 0.9942327
#>
#> $wRMSD
#> [1] 0.01077901
```

Bending the same correlation matrix using the unweighted Schaeffer (2014):

```r
bend(V2, method="lrs")
#> Unweighted bending
#> max.iter = 10000
#> method = lrs
#> Found a correlation matrix.
#> Convergence met after 39 iterations.
#> $bent
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.8935779 0.7794136 0.4716892 0.3599505
#> [2,] 0.8935779 1.0000000 0.9128285 0.7211811 0.4716892
#> [3,] 0.7794136 0.9128285 1.0000000 0.9128285 0.7794136
#> [4,] 0.4716892 0.7211811 0.9128285 1.0000000 0.8935779
#> [5,] 0.3599505 0.4716892 0.7794136 0.8935779 1.0000000
#>
#> $init.ev
#> [1]  3.99475997  0.98523500  0.23646897 -0.03122893 -0.18523500
#>
#> $final.ev
#> [1] 3.9083638707 0.9183589993 0.1656638620 0.0071038704 0.0005093975
#>
#> $min.dev
#> [1] -0.07881887
#>
#> $max.dev
#> [1] 0.07168916
#>
#> $loc.min.dev
#> row col
#>   4   2
#>
#> $loc.max.dev
#> row col
#>   5   2
#>
#> $ave.dev
#> [1] -0.02038503
#>
#> $AAD
#> [1] 0.04906069
#>
#> $Cor
#> [1] 0.9854043
#>
#> $RMSD
#> [1] 0.05298397
```

---

## References

Jorjani, H., Klie. L., & Emanuelson, U. (2000). A simple method for weighted bending of genetic (co)variance matrices. _J. Dairy Sci._ 86(2): 677--679. [doi:10.3168/jds.S0022-0302(03)73646-7](<https://doi.org/10.3168/jds.S0022-0302(03)73646-7>)

Schaeffer, L. R. (2014). Making covariance matrices positive definite. Available at: <http://animalbiosciences.uoguelph.ca/~lrs/ELARES/PDforce.pdf>
