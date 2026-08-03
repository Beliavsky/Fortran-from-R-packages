# R4GoodPersonalFinances (development version)

# R4GoodPersonalFinances 1.2.0

## Minor changes

* Changed implementation of formula for expected utility in the `calc_expected_utility()` internal function. Now the Constant Relative Risk Aversion (CRRA) utility and the second derivative of the CRRA utility are implemented in a way that provide more stable gradient for calculating optimal portfolio allocation.
As a result of this change, the optimal allocation no longer jumps unexpectedly and is more stable over time.
* Added new and improved existing debug messages visible when `debug = TRUE` argument is passed to `simulate_scenario()` function.
* The function `calc_gompertz_mode()` is now exported and it can be used to calculate the mode of the Gompertz distribution based on expected `life_expectancy`.

## Bug fixes

* Fixed discretionary spending dropping unexpectedly to zero in some edge cases caused by division by zero while calculating ratio in the `calc_discretionary_spending()` function. Now, if the ratio is undefined, it it set to 1.
* Fixed undefined Gompertz survival probability in a case when `current_age` equals `target_age` and `max_age`. Now, in such cases, the survival probability is set explicitly to 0.
* Fixed discretionary spending at the last year of a simulation. Now the discretionary spending is set to net-worth at the last year of a simulation which results in zero of financial wealth at the end of the last year.


# R4GoodPersonalFinances 1.1.0

This version improve the performance of portfolio optimization algorithms,
ability to parallelize the Monte Carlo simulations, and simulation caching.
It brings also new useful functions as `render_scenario_snapshot()`,  `plot_future_saving_rates()`, and `get_default_gompertz_parameters()`.
With `render_scenario_snapshot()` function, it is possible to get a table with cash flow and net worth statements for a given year in a simulated scenario.
The result of new features in this version can be seen in the blog post: [Can Juliet manage without Romeo? How much she should spend, save, and invest to achieve that?](https://www.r4good.academy/en/blog/financial-report-individual-minimal-example/index.en.html).


## Minor changes

* Changed default portfolio optimization algorithm to `NLOPT_LD_SLSQP`
(Sequential Least-Squares Quadratic Programming). The optimization is now more robust and faster. 
* Added ability to pass `opts` list to `simulate_scenario()` function with optimization parameters `nloptr` optimization function and thus changing algorithm or its options. The list of possible options and their default values can be obtain with `nloptr::nloptr.print.options()`. 
* Added ability to pass a single seed value or vector of seeds `simulate_scenario()` function to reproduce the same results for the same scenario, also for Monte Carlo samples.
* Added `auto_parallel` argument to `simulate_scenario()` function to automatically parallelize the Monte Carlo simulations.
* Changed behavior of `use_cache` argument in `simulate_scenarios()` function 
to caching also all data of all Monte Carlo samples for a scenario (not only each sample separately). 
* Changed the reordering behavior in the `plot_scenarios()` which does not reorder scenarios anymore.
* Changed optimal allocation for the last year of a simulation. 
It is now overwritten with the optimal allocation from the previous year,
to avoid the edge case when net-worth is going to zero.
* Added `render_scenario_snapshot()` function to render a table with cash flow and net worth statements for a given year (first by default) in a simulated scenario.
* Added `plot_future_saving_rates()` for plotting the future saving rates.
* Added print methods for `Household` and `HouseholdMember` objects.
* Added ability to plot expected allocations based on the aggregated results from Monte Carlo samples with `plot_expected_allocation()`.
* Added `age()` function as convenient wrapper for getting age of a household member when defining triggers for expected income or spending streams.
* Added `get_default_gompertz_parameters()` function to get default parameters for the Gompertz model for a given country, sex, and age.
* Added, changed, or fixed multiple minor details in plots.
* Updated default capital market assumptions.


# R4GoodPersonalFinances 1.0.0

This major release introduces probably first in the world (open-sourced) implementation of a new multilevel life-cycle modeling of a household finances. 
It connects life-cycle models with single-period net-worth Markowitz 
optimization models and allows for high degree of personalization.
The implementation is based on a novel theoretical framework, 
described by Thomas M. Idzorek and Paul D. Kaplan in 
"Lifetime Financial Advice: A Personalized Optimal Multilevel Approach" (2024).

## Major changes

* Added `Household` R6 class that aggregates information about a household 
and its members and  `HouseholdMember` class that aggregates information 
about a single member of a household.
* Added functions responsible for setting, resetting and getting information about disk cache for caching the results of the simulations: `set_cache()`, `reset_cache()` and `get_cache()`.
* Added `calc_life_expectancy()` function to calculate the life expectancy of a person based on the Gompertz model.
* Added `plot_life_expectancy()` function to visualize the probability of dying and life expectancy of every member of a `Household` object and similar `plot_survival()` function for visualizing the probability of survival of each of the members of a `Household` and the probability of at least one member of a `Household` to survive.
* Added `create_portfolio_template()` function to create an exemplary `Portoflio` object - a template for default portfolio with two asset classes: `GlobalStocksIndexFund` and `InflationProtectedBonds` and with some reasonable default values that should be updated and customised further.
* Added `calc_portfolio_parameters()` function to calculate the parameters of a `Portfolio` object like current value, weights of assets, expected return and standard deviation of the whole portfolio.
* Added `calc_optimal_asset_allocation()` function to calculate optimal allocation for a particular `Household` and `Portfolio` objects.
* Added `plot_optimal_portfolio()` function to visualize the current and optimal allocation to assets for a particular `Portfolio` object.
* Added `simulate_scenario()` function to simulate a single scenario of cash flows for a particular `Household` and `Portfolio` objects using expected returns of the assets in portfolio or random returns in multiple Monte Carlo samples.
* Added `simulate_scenarios()` function to simulate multiple scenarios of cash flows for a particular `Household` and `Portfolio` objects. The scenarios parameters can triggered by different events like different years of retirement.
* Added `plot_scenarios()` function to visualize the scenarios metrics that help with choosing the optimal parameters for the best scenario, e.g. the optimal year for retirement.
* Added `plot_expected_allocation()` function to visualize the expected allocation of assets for a particular scenario over the household life cycle.
* Added `plot_expected_capital()` function to visualize the expected financial, human, and total capital, as well as liabilities for a particular scenario over the household life cycle.
* Added `plot_future_income()` and `plot_future_spending()` functions to visualize the future income and spending for a particular scenario over the household life cycle. The `plot_future_spending()` function allows to plot the relation of discretionary and non-discretionary spending, structure of non-discretionary spending and possible levels of discretionary spending over time based on Monte Carlo simulations.


# R4GoodPersonalFinances 0.2.0

This release provides functions for calculating probability of retirement ruin and functions for calibrating the underlying Gompertz model with the use of the mortality rates.

## Major changes

* Added `calc_retirement_ruin()` and `plot_retirement_ruin()` to calculate and visualize the probability of retirement ruin.
* Added `retirement-ruin` option to `run_app()` function for
running package new interactive app.
* Added `calc_gompertz_survival_probability()` to calculate the probability of survival using the Gompertz model. 
* Added `calc_gompertz_parameters()` to calculate the parameters of the Gompertz model based on mortality rates and `calc_gompertz_joint_parameters()` to calculate the parameters of the joint Gompertz model.
* Added `plot_gompertz_calibration()` for visually checking the calibration of the Gompertz model and `plot_joint_survival()` for plotting the survival probabilities for two individuals and their joint survival probability.
* Added `read_hmd_life_tables()` for reading the 'Human Mortality Database' (HMD) life tables files.

# R4GoodPersonalFinances 0.1.1

## Bugfixes

* Fixed apps errors while executing locally.

# R4GoodPersonalFinances 0.1.0

Initial release of the package includes functions designed to 
support the optimal allocation of your investment portfolio
with the use of the 'Merton Share' formula.

* Added `calc_purchasing_power()` and `plot_purchasing_power()` to calculate and visualize the effect of real interest rates—whether positive (compounding) or negative (inflation)—on the purchasing power of savings.  
* Added `calc_optimal_risky_asset_allocation()` for calculating optimal allocation using the 'Merton Share' formula.  
* Added `calc_risk_adjusted_return()` and `plot_risk_adjusted_returns()` to calculate and visualize the effect of allocation to risky assets on the risk-adjusted return.
* Added `run_app()` function that runs package interactive apps:
`risk-adjusted-returns` and `purchasing-power`.
