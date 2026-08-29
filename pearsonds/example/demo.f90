program demo
use pearsonds_mod
implicit none
type(pearson_params_t) :: pars, recovered
real(kind=dp) :: moments(4)
real(kind=dp), allocatable :: x(:), dens(:), cdf(:)

pars%family = pearson_type_iv
pars%npar = 4
pars%par = [6.0_dp, 2.0_dp, 0.5_dp, 1.25_dp]

x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
dens = dpearson(x, pars)
cdf = ppearson(x, pars)
moments = pearson_moments(pars)
recovered = pearson_fit_m(moments(1), moments(2), moments(3), moments(4))

write(*,'(a,*(1x,g0))') 'density:', dens
write(*,'(a,*(1x,g0))') 'cdf:    ', cdf
write(*,'(a,*(1x,g0))') 'moments:', moments
write(*,'(a,i0)') 'recovered family: ', recovered%family
end program demo
