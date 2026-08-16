program test_core
   use vgam
   implicit none
   integer :: failures
   failures = 0

   call check_close("normal cdf", pnorm_v(1.2_dp,0.0_dp,1.0_dp), &
      0.8849303297782918_dp, 2.0e-11_dp, failures)
   call check_close("normal quantile", qnorm_v(0.975_dp,0.0_dp,1.0_dp), &
      1.959963984540054_dp, 5.0e-9_dp, failures)
   call check_close("gamma cdf", pgamma_v(2.3_dp,1.7_dp,0.8_dp), &
      0.8407159904758339_dp, 2.0e-10_dp, failures)
   call check_close("beta cdf", pbeta_v(0.35_dp,2.4_dp,3.2_dp), &
      0.3722253423526247_dp, 2.0e-10_dp, failures)
   call check_close("poisson cdf", ppois_v(5,3.4_dp), &
      0.8705423698828116_dp, 2.0e-10_dp, failures)
   call check_close("binomial cdf", pbinom_v(3,10,0.2_dp), &
      0.8791261184_dp, 2.0e-12_dp, failures)
   call check_close("negative binomial pmf", dnbinom_v(4,2.5_dp,0.4_dp), &
      0.1183387545988211_dp, 2.0e-12_dp, failures)
   call check_close("logit inverse", link_inverse(link_value(0.37_dp,link_logit),link_logit), &
      0.37_dp, 2.0e-14_dp, failures)
   call check_close("probit inverse", link_inverse(link_value(0.37_dp,link_probit),link_probit), &
      0.37_dp, 2.0e-9_dp, failures)
   call check_close("gumbel inversion", pgumbel(qgumbel(0.23_dp,1.2_dp,0.7_dp),1.2_dp,0.7_dp), &
      0.23_dp, 2.0e-12_dp, failures)
   call check_close("rayleigh inversion", prayleigh(qrayleigh(0.73_dp,1.4_dp),1.4_dp), &
      0.73_dp, 2.0e-12_dp, failures)
   call check_close("pareto4 inversion", ppareto4(qpareto4(0.61_dp,1.5_dp,2.1_dp,0.8_dp,0.3_dp), &
      1.5_dp,2.1_dp,0.8_dp,0.3_dp), 0.61_dp, 2.0e-11_dp, failures)
   call check_close("inverse Gaussian inversion", &
      pinvgaussian(qinvgaussian(0.42_dp,1.3_dp,2.2_dp),1.3_dp,2.2_dp), &
      0.42_dp, 2.0e-8_dp, failures)
   call check_close("Gompertz inversion", pgompertz(qgompertz(0.7_dp,0.15_dp,0.8_dp),0.15_dp,0.8_dp), &
      0.7_dp, 2.0e-10_dp, failures)
   call check_close("Lindley inversion", plindley(qlindley(0.4_dp,1.3_dp),1.3_dp), &
      0.4_dp, 2.0e-9_dp, failures)
   call check_close("Lambert W", lambert_w0(1.0_dp)*exp(lambert_w0(1.0_dp)), &
      1.0_dp, 2.0e-12_dp, failures)
   call check_close("GEV inversion", pgev(qgev(0.36_dp,0.4_dp,1.2_dp,0.15_dp),0.4_dp,1.2_dp,0.15_dp), &
      0.36_dp, 2.0e-12_dp, failures)
   call check_close("GPD inversion", pgpd(qgpd(0.66_dp,0.2_dp,0.9_dp,-0.1_dp),0.2_dp,0.9_dp,-0.1_dp), &
      0.66_dp, 2.0e-12_dp, failures)
   call check_close("Kumaraswamy inversion", pkumar(qkumar(0.58_dp,1.7_dp,2.3_dp),1.7_dp,2.3_dp), &
      0.58_dp, 2.0e-12_dp, failures)

   if (failures /= 0) then
      print '(a,i0)', "test_core failures: ", failures
      error stop 1
   end if
   print '(a)', "test_core: PASS"

contains
   subroutine check_close(name, actual, expected, tol, failures)
      character(*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tol
      integer, intent(inout) :: failures
      if (abs(actual-expected) > tol .or. .not.(abs(actual) <= huge(actual))) then
         print '(a,2(1x,es24.16),a,es12.4)', trim(name)//" FAIL",actual,expected," tol=",tol
         failures=failures+1
      end if
   end subroutine check_close
end program test_core
