program wolf_die
    use pmultinom_module, only : dp, pmultinom
    implicit none

    real(dp) :: p_fair, p_model
    real(dp) :: fair(6), theoretical(6), upper(6)

    fair = 1.0_dp / 6.0_dp
    theoretical = [0.17649_dp, 0.17542_dp, 0.15276_dp, 0.15184_dp, 0.17227_dp, 0.17122_dp]
    upper = 3630.0_dp

    p_fair = 1.0_dp - pmultinom(20000, fair, upper=upper)
    p_model = 1.0_dp - pmultinom(20000, theoretical, upper=upper)

    print '(a,es14.6)', 'fair-die tail probability: ', p_fair
    print '(a,es14.6)', 'Fougere-model tail probability: ', p_model
end program wolf_die
