program basic
    use actuar
    implicit none
    type(aggregate_dist_t) :: aggregate
    real(dp) :: severity(5)

    severity = [0.0_dp, 0.45_dp, 0.30_dp, 0.15_dp, 0.10_dp]
    aggregate = panjer_poisson(severity, 3.0_dp, tol=1.0e-10_dp)
    aggregate%x_scale = 1000.0_dp

    print '(a,f10.4)', 'Pareto 99% quantile: ', qpareto(0.99_dp,3.0_dp,1000.0_dp)
    print '(a,f10.4)', 'Aggregate mean:       ', aggregate%mean()
    print '(a,f10.4)', 'Aggregate 95% VaR:    ', aggregate_var(aggregate,0.95_dp)
    print '(a,f10.4)', 'Aggregate 95% CTE:    ', aggregate_cte(aggregate,0.95_dp)
end program basic
