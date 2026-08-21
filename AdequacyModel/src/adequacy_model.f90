! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_model
    use adequacy_kinds, only: dp
    use adequacy_interfaces, only: objective_fn, density_fn, cdf_fn
    use adequacy_stats, only: descriptive_result, descriptive, ttt_curve
    use adequacy_optim, only: optimize_result, pso_optimize, nelder_mead_optimize
    use adequacy_optim, only: bfgs_optimize, cg_optimize, sann_optimize
    use adequacy_gof, only: goodness_result, goodness_fit, goodness_from_mle
    implicit none
    public
end module adequacy_model
