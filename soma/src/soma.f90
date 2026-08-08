! SPDX-License-Identifier: GPL-2.0-only
module soma
    use soma_kinds, only : dp
    use soma_types, only : soma_bounds, soma_options, soma_result, soma_cost_function, &
                           make_bounds, bounds, all2one, t3a, pareto, &
                           strategy_all2one, strategy_t3a, strategy_pareto
    use soma_random, only : soma_set_seed
    use soma_optimizer, only : soma_optimize
    implicit none
    public
end module soma
