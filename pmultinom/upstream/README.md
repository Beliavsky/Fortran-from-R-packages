`pmultinom` is a library for calculating multinomial probabilities. The probabilities that can be calculated include the multinomial cumulative distribution function:
$$P(N_1 \le u_1, N_2 \le u_2, \cdots, N_k \le u_k)$$
In this case the usage would be
```
pmultinom(upper=us, size=n, probs=ps, method="exact")
```
where `us` is the vector containing $u_1, u_2, \cdots, u_k$, and `n` and `ps` are the parameters of the multinomial distribution. This usage is analogous to the use of `pbinom`. Another important case is the probability of seeing more than some minimum number of observations in each category:
$$P(N_1 > l_1, N_2 > l_2, \cdots, N_k > l_k)$$
In this case the usage would be
```
pmultinom(lower=ls, size=n, probs=ps, method="exact")
```
where this time `ls` is the vector containing $l_1, l_2, \cdots, l_k$. Notice that in this case these are greater than signs, not greater than or equal signs. This is analogous to the usage of `pbinom` with `lower.tail=FALSE`. With some creativity, these can be adapted to calculate the probability that the maximum or minimum of a multinomial random vector is a given number, or that a given category will be the most or least observed. `pmultinom` also supports a more general usage, in which both lower and upper bounds are specified:
$$P(l_1 < N_1 \le u_1, l_2 < N_2 \le u_2, \cdots, l_k < N_k \le u_k)$$
In this case the usage would be
```
pmultinom(lower=ls, upper=us, size=n, probs=ps, method="exact")
```

See `vignette("pmultinom")` for the above text in Latex plus an example application. Many thanks to Aislyn Schalck for advice and encouragement.
