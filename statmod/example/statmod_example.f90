program statmod_example
use statmod
use r_compat, only: dp
implicit none
type(quad_rule_t) :: g
real(dp) :: m, v
print '(a,es16.8)', 'P(IG <= 1; mean=2, dispersion=.5) = ', &
   pinvgauss(1.0_dp, mean=2.0_dp, dispersion=0.5_dp)
g=gauss_quad(5,'legendre')
print '(a,5f12.7)', 'Legendre nodes: ',g%nodes
call expected_deviance(2.0_dp,'poisson',m=m,v=v)
print '(a,2f14.9)', 'Poisson expected deviance mean/variance: ',m,v
end program
