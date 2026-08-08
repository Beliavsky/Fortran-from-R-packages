# optimflex 0.1.8
* Added `levenberg_marquardt()`: a Levenberg-Marquardt optimizer 
* Fixed a bug in which `lower` and `upper` bounds were accepted but not
  applied during optimization.
* Fixed a bug in which some convergence criteria were accepted but not
  evaluated.

# optimflex 0.1.7
* Fixed bugs in dogleg and double_dogleg so that when the user supplies an exact Hessian, it is used directly in the optimization rather than being overridden by the BFGS approximation.

# optimflex 0.1.6
* Replaced the quadratic example with the Rosenbrock example. Accordingly, removed dfp.R, which had failed to converge on this 
example in certain environments.

# optimflex 0.1.5

Bug Fixes

* Fixed a convergence detection failure when the optimizer reaches the exact
solution in a single step. A post-step gradient check has been added to all
eight optimization algorithms.
* Replaced test-rosenbrock.R with test-quadratic.R to verify convergence
on a simple quadratic function, consistent with the example in the dfp
documentation.

# optimflex 0.1.4

* Resubmission for CRAN.
* Removed tests/testthat/test-rosenbrock.R, which was inadvertently included in the previous submission and caused a convergence failure on M1 Mac (ATLAS BLAS)

# optimflex 0.1.3

* Resubmission for CRAN.
* Replaced the example in 'dfp.R' with a simpler one.

# optimflex 0.1.1

* Resubmission for CRAN.
* Fixed CRAN check errors in `multirel` examples and refined LaTeX documentation in `dfp.R`.
* Implements a flexible suite of optimization algorithms (BFGS, DFP, Dogleg, Newton, etc.).
* Features strict convergence control with eight combinable criteria.
* Includes a built-in mechanism for final Hessian positive definiteness verification.