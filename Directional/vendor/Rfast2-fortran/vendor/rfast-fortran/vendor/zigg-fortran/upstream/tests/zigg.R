
library(zigg)
options(digits=6, scipen=99)

N <- 1e5
zsetseed(12345)

rn <- zrnorm(N)
round(summary(rn), digits=4)

re <- zrexp(N)
round(summary(re), digits=4)

ru <- zrunif(N)
round(summary(ru), digits=4)

zsetseed(12345)
all.equal(rn, zrnorm(N))
all.equal(re, zrexp(N))
all.equal(ru, zrunif(N))
