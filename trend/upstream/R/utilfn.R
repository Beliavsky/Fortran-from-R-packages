##   file utilfn.R part of package trend
##   Copyright (C) 2015-2018 Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
# Internal functions

# calculates varS for MK-test and Sen's slope
.varmk <- function(t, n)
{
    tadjs <- sum(t * (t - 1) * (2 * t + 5))
    varS <- (n * (n-1) * (2 * n + 5) - tadjs) / 18
    return(varS)
}

## calculates D for tau estimation
.Dfn <- function(t, n)
{
    tadjd <- sum(t * (t - 1))
    D <- sqrt(1/2 * n * (n - 1) - 1/2 * tadjd) *
        sqrt(1/2 * n * (n - 1))
    return(D)
}

## caclulates the Mann-Kendall Score
## .mkScore <- function(x){
##    n <- length(x)
##    S <- 0.0
##  
##    for (k in 1:(n-1))
##    {
##        for (j in k:n)
##        {
##            S <- S + sign(x[j] - x[k]);
##        }
##    }
##    return(S)
##}

## recoded as recommended from
## Martin Schmelzer, email from 2017-11-11
.mkScore <- function(x){
    n <- length(x)
    S <- 0.0   
    for(j in 1:n) {
        S <- S + sum(sign(x[j] - x[1:j]))
    }	
    return(S)
}



## Calculate K for partial Mann-Kendall Trend Test
.K <- function(x, z)
{
    n <- length(x)
    K <- 0.0

    for (i in 1:(n-1)){
        for (j in i:n){
            K <- K + sign((x[j] - x[i]) * (z[j] - z[i]))
        }
    }
    return(K)
}

## Calculate R
.R <- function(x)
{
    n <- length(x)
    R <- numeric(n)
    for (j in 1:n)
    {
        s <- 0.0
        for (i in 1:n)
        {
            s <- s + sign(x[j] - x[i])
        }

        R[j] <- (n + 1 + s ) / 2
    }
    return(R)
}
