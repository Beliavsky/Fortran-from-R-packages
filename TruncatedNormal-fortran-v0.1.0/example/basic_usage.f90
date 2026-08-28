program basic_usage
 use truncated_normal
 use r_mod, only: set_seed_int
 implicit none
 real(dp) :: mu(2), sigma(2,2), lower(2), upper(2)
 real(dp), allocatable :: draws(:,:)
 type(prob_result) :: p

 mu = [0.0_dp, 0.0_dp]
 sigma = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2,2])
 lower = [-1.0_dp, -1.0_dp]
 upper = [1.0_dp, 1.0_dp]

 call set_seed_int(1234)
 p = pmvnorm(mu, sigma, lower, upper, B=20000, qmc=.true.)
 print '(a,f12.8,a,es12.4)', 'P(lower < X < upper) = ', p%prob, '  relative error = ', p%relerr

 draws = rtmvnorm(5, mu, sigma, lower, upper)
 print '(a)', 'Five truncated-normal draws:'
 print '(2f14.7)', transpose(draws)
end program basic_usage
