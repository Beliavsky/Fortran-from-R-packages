program test_fit
   use skewunit
   implicit none
   type(skewunit_fit_result) :: fit
   type(skewunit_choice_result) :: choice
   real(dp), parameter :: x1(35) = [ &
      0.16518712992191037_dp,0.20660428701775424_dp,0.69030478958852526_dp, &
      0.26859060902701748_dp,0.11827447856213233_dp,0.6328254894362334_dp, &
      6.383051949928375e-05_dp,0.012037823644668795_dp,0.32279263524283247_dp, &
      0.17429586178682824_dp,0.42386706695502419_dp,0.22561623630905991_dp, &
      0.041826091258883119_dp,0.60132159371355087_dp,0.90285028305098824_dp, &
      0.0084685662636990854_dp,0.79595008086626917_dp,0.0021334623477326895_dp, &
      0.019874007109502528_dp,0.086579758336783094_dp,0.58175616571661304_dp, &
      0.030446611440016702_dp,0.20424978862350363_dp,0.31100488411965171_dp, &
      0.18826993600235953_dp,0.14300407187503184_dp,0.024601959787036356_dp, &
      0.75285564758729573_dp,0.99999603641359092_dp,0.76804608705308897_dp, &
      0.61126138999537516_dp,0.30636698158029507_dp,0.11675989483583479_dp, &
      0.59348909359560031_dp,0.038183270531717141_dp ]
   real(dp), parameter :: x2(20) = [ &
      0.77447775847187006_dp,0.098516776276084267_dp,0.70476008986547367_dp, &
      0.44518613727883577_dp,0.59548422583169347_dp,0.14730797522827838_dp, &
      0.31595549404613554_dp,0.21232839064024009_dp,0.080607221323583569_dp, &
      0.71237294071678658_dp,0.81689044183578008_dp,0.70165342598473079_dp, &
      0.8185039544283591_dp,0.38745573014175733_dp,0.053966965665306052_dp, &
      0.057432435613069883_dp,0.24553721333962636_dp,0.28063855909242175_dp, &
      0.60224364736629732_dp,0.60228405268687601_dp ]
   integer :: i

   call estimate_skewunit(x1,family_asin,family_jsb,fit,est_var=.true.)
   if (fit%convergence /= 0) error stop 'case 2 fit did not converge'
   call check_close('lambda fit',fit%coefficients(1),-0.8085452505617611_dp,3.0e-4_dp)
   call check_close('delta fit',fit%coefficients(2),0.5399003329126715_dp,3.0e-4_dp)
   call check_close('logLik fit',fit%loglik,14.748761567261829_dp,2.0e-6_dp)
   if (.not.fit%std_error_available) error stop 'standard errors unavailable'
   if (any(fit%std_error(1:2) <= 0.0_dp)) error stop 'invalid standard errors'

   call estimate_skewunit(x2,family_sbeta,family_none,fit,est_var=.false.)
   if (fit%convergence /= 0) error stop 'baseline shape fit did not converge'
   if (fit%coefficients(1) <= 0.0_dp) error stop 'baseline delta nonpositive'

   call choose_skewunit(x2,choice,criteria='AIC',est_var=.false.,maxit=3000)
   do i = 2, 30
      if (choice%summary(i)%criterion < choice%summary(i-1)%criterion) &
         error stop 'model summary not sorted'
   end do
   if (abs(choice%best_fit%aic-choice%summary(1)%criterion) > 1.0e-7_dp) &
      error stop 'best model mismatch'

   print '(a)', 'test_fit: PASS'

contains

   subroutine check_close(name, got, expected, tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got, expected, tol
      if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
         print '(a,2(1x,es24.16))', trim(name)//' got/expected:',got,expected
         error stop 1
      end if
   end subroutine check_close

end program test_fit
