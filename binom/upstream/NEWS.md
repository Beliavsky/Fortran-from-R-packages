# binom 1.1-2

* Fixed `binom.confint()` with `method = "prop.test"` ignoring the
  `conf.level` argument: `stats::prop.test()` is now called with
  `conf.level = conf.level`, so intervals are no longer always
  computed at the default 95% level.
* Arianna D. Cascone is now the maintainer; Sundar Dorai-Raj remains
  an author.

# binom 1.1-1.1

* CRAN maintenance release (2022-05-02).

# binom 1.1-1

* CRAN maintenance release (2014-01-01).

# binom 1.0-5

* CRAN maintenance release (2009-05-26).

# binom 1.0-4

* CRAN maintenance release (2008-10-03).

# binom 1.0-3

* CRAN maintenance release (2008-01-02).

# binom 1.0-2

* CRAN maintenance release (2007-11-17).

# binom 1.0-1

* CRAN maintenance release (2007-02-13).

# binom 1.0-0

* Initial CRAN release (2006-08-09).
