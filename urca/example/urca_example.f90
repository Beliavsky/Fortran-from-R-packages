program urca_example
   use urca, only : dp, ur_test_result, johansen_result, adf_test, johansen_test, &
      UR_DRIFT, LAG_AIC, JO_TRACE, JO_CONST, JO_LONGRUN
   implicit none
   integer, parameter :: n = 300
   real(dp) :: y(n), x(n,2), e1, e2, state
   integer :: i, nseed
   integer, allocatable :: seed(:)
   type(ur_test_result) :: adf
   type(johansen_result) :: jo

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed=20260828
   call random_seed(put=seed)

   y(1)=0.0_dp
   x(1,:)=0.0_dp
   state=0.0_dp
   do i=2,n
      call normal_random(e1)
      call normal_random(e2)
      y(i)=0.65_dp*y(i-1)+e1
      state=0.5_dp*state+0.5_dp*e2
      x(i,1)=x(i-1,1)+e1
      x(i,2)=x(i,1)+state
   end do

   adf=adf_test(y,UR_DRIFT,6,LAG_AIC)
   jo=johansen_test(x,JO_TRACE,JO_CONST,2,JO_LONGRUN)

   print '(a,f10.4,a,i0)', 'ADF tau = ',adf%statistic(1),', selected lag = ',adf%lags
   print '(a,*(f10.4,1x))', 'Johansen eigenvalues = ',jo%lambda
   print '(a,*(f10.4,1x))', 'Johansen trace statistics = ',jo%teststat
contains
   subroutine normal_random(z)
      real(dp), intent(out) :: z
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end subroutine normal_random
end program urca_example
