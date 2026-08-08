# Validation

The automated tests exercise:

- Hooke-Jeeves convergence on the Rosenbrock function
- bounded Hooke-Jeeves with an active upper bound
- modified Nelder-Mead on Rosenbrock
- transformed-bound Nelder-Mead with an active bound
- lower-only and upper-only transformed variables
- MADS on a nonsmooth absolute-value objective
- MADS full polling and iteration-log output
- maximization
- evaluation-limit and invalid-input statuses
- polymorphic user-data callbacks
- monitor cancellation

Builds are tested in two modes:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

and

```text
-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface
```

The package test functions use independently known solutions, including the
Rosenbrock minimum `(1,1)`, active-bound quadratic minima, and a nonsmooth
absolute-value minimum.
