# ecd 0.9.2.3

2022-05-10: Removed RcppFaddeeva dependency.

2022-05-02: Removed fOption dependency. The research from 2017-2022 revealed that the core structure of the lambda distribution is a volatility kernel based on the Stable Count Distribution, which is documented formally in Wikipedia:
https://en.wikipedia.org/wiki/Stable_count_distribution

# ecd 0.8.3

2017-01-05: Primarily a fix for RSQLite. Functions supporting VIX option pricing model are added.
The paper on PDE and Local Volatility Function is about to be finished.

# ecd 0.8.2

2016-07-09: Finalized the quartic lambda model and the supporting library for the paper.
The paper "Closed Form Solution and Term Structure for SPX Options"
is available at http://ssrn.com/abstract=2805769
 
# ecd 0.7.0

2016-04-30: The quartic lambda option pricing model (lambda=4) is added. 
It provides computation engine for volatility smile of SPX options. 
The distribution now can take a real number for lambda (instead of a rational number).
This allows user to use optimx to fit the market data to the pricing model.

# ecd 0.6.4

2015-12-23: This is the first official release that passed the CRAN policy.
The extdata folder has been trimmed down to fit the 5MB limit.
The heavy computing test cases were moved out so that check can finish within 20 minutes.
The two companion papers on the elliptic distribution and elliptic option pricing model
were published to SSRN as well.

# ecd 0.6.0

This release is ready for CRAN on 2015-12-01. 
It includes all the functions for ecd and ecld.

# ecd 0.1.0 

The first ‘ecd’ release was made on 2015-05-31




