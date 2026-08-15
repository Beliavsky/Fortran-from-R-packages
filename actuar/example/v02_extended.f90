program v02_extended
    use actuar
    implicit none
    real(dp) :: sev(3), freq(5)
    real(dp) :: pc(1), tc(1,1), pw(1), tw(1,1)
    type(aggregate_dist_t) :: agg
    type(ruin_result_t) :: ruin

    print '(a,f12.6)', 'Feller-Pareto P(X <= 6): ', &
        pfpareto(6.0_dp,1.0_dp,2.5_dp,1.4_dp,1.7_dp,3.0_dp)

    sev = [0.15_dp,0.55_dp,0.30_dp]
    freq = [0.10_dp,0.20_dp,0.30_dp,0.25_dp,0.15_dp]
    agg = aggregate_exact(sev,freq,1000.0_dp)
    print '(a,f12.3)', 'Exact aggregate mean: ', agg%mean()
    print '(a,f12.3)', 'Exact aggregate 95% VaR: ', aggregate_var(agg,0.95_dp)

    pc = 1.0_dp
    tc(1,1) = -2.0_dp
    pw = 1.0_dp
    tw(1,1) = -0.5_dp
    ruin = ruin_phase_type(pc,tc,pw,tw,1.0_dp)
    print '(a,f12.6)', 'Ruin probability at surplus 1: ', ruin%probability(1.0_dp)
end program v02_extended
