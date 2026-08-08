# neighbours

Neighbourhood functions are key components of
local-search algorithms such as Simulated Annealing or
Threshold Accepting.  These functions take a solution
and return a slightly-modified copy of it, i.e. a
neighbour.  The package provides a function,
neighbourfun(), that produce a neighbourhood function,
based on various parameters such as minimum or maximum
values.  Currently supported are numeric and logical
solutions.  The algorithms were originally created for
portfolio-optimization models, but can also be used for
other models.

[ [More] ](https://enricoschumann.net/R/packages/neighbours/)

## Installing the package

The latest released version is available from
http://enricoschumann.net. In an R session, just type:

    install.packages('neighbours',
                     repos = c('https://enricoschumann.net/R',
                               getOption('repos')))
