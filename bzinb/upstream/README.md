# BZINB
Bivariate Zero-Inflated Negative Binomial Model Estimation  

This package is based on a draft paper ``A Bivariate Zero-Inflated Negative Binomial Model For Identifying Underlying Dependence" (Hunyong Cho, John Preisser, Chuwen Liu \& Di Wu 2019+ (In preparation)).  

See the following toy example for fun.

```{r}
library(bzinb)

# generating n x 2 matrix (two vectors)
set.seed(2)
data1 <- rbzinb(n = 20, a0 = 1, a1 = 2, a2 = 1,
               b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2,
               p3 = 0.2, p4 = 0.1)

# getting the underlying correlation (rho) through maximum likelihood estimate.
bzinb(xvec = data1[,1], yvec = data1[,2], showFlag = F)


# generating (additional two vectors)
set.seed(3)
data2 <- rbzinb(n = 20, a0 = 2, a1 = 1, a2 = 1, 
                b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2, 
                p3 = 0.2, p4 = 0.1)
data3 <- t(cbind(data1, data2))
pairwise.bzinb(data3, showFlag = TRUE)

```
