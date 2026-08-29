library(DPQ)
### Testing the Fortran / C code

###  originally == /u/maechler/R/MM/NUMERICS/dpq-functions/wienergerm-pchisq-tst.R

Aeq <- function(x,y) all.equal(x,y, tol = 1e-10)

for (Ftn in c(TRUE,FALSE)) {
    cat("Fortran =",Ftn, "\n")
    stopifnot(
              Aeq(pchisqW.(print(7), df=1, ncp=2, Fort=Ftn)[[1]]$p,
                  0.890202842),
              Aeq(pchisqW.(7, df=1, ncp=2, Fort=Ftn, vari=print('f'))[[1]]$p,
                  0.8908801586),
        ##
              Aeq(pchisqW.(print(44), df=4, ncp=16, Fort=Ftn)[[1]]$p,
                  0.99045176549),
        ##
              Aeq(pchisqW.(44, df=4, ncp=16, Fort=Ftn, vari=print('f'))[[1]]$p,
                  0.9904601514265),
        ##
              Aeq(pchisqW.(print(50), df=3.5, ncp=20, Fort=Ftn)[[1]]$p,
                  0.9912862375384),
              Aeq(pchisqW.(50, df=3.5,ncp=20, Fort=Ftn, vari=print('f'))[[1]]$p,
                  0.9912928203989)
              )
}

pchisqW(1:3, df=1, ncp=2)

x1 <- c(1:3,10,50,100,199:205) # problem at  x == 201 = df+ncp  !!!
p1 <- pchisqW(x1, df=200, ncp=1)
Aeq(p1,
    c(3.132799411e-189, 2.4267526025e-159, 6.029138178e-142, 3.7252821178e-91,
      8.4350675182e-30, 2.4797380646e-10, 0.47350428951, 0.49343195544,
      1, ## was NaN
      0.53307133169, 0.55270012044, 0.57213676726, 0.59133736744))
## [1] "Mean relative difference: 3.384587e-10"
