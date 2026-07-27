# tvm 0.4

* Add year fraction and compounding arguments to both `xnpv` and `xirr`, providing greater flexibility

# tvm 0.3

* Add vignette and NEWS
* Prevent negative discount factors when calculating them from swap curves
* Split zero rates in zero nominal and zero effective
* A functor argument is added to the `rate_curve` constructor, to allow the user to specify how the interpolation should be performed.
* Changing the default spline interpolation method from `natural` to `monoH.FC`, as to respect monotonicity if present
* Add new irregular functions, `xnpv` and `xirr`
