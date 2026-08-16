program test_moments
   use discrete_weibull
   implicit none
   real(dp) :: q,c,m1,m2,v
   integer :: fails

   fails = 0
   q = 0.4_dp
   m1 = Edweibull(q,1.0_dp,zero=.false.)
   m2 = E2dweibull(q,1.0_dp,zero=.false.)
   v = Vdweibull(q,1.0_dp,zero=.false.)
   if (abs(m1-1.0_dp/(1.0_dp-q)) > 1.0e-14_dp) fails=fails+1
   if (abs(m2-(1.0_dp+q)/(1.0_dp-q)**2) > 1.0e-14_dp) fails=fails+1
   if (abs(v-q/(1.0_dp-q)**2) > 1.0e-14_dp) fails=fails+1

   m1 = Edweibull(q,1.0_dp,zero=.true.)
   if (abs(m1-q/(1.0_dp-q)) > 1.0e-14_dp) fails=fails+1

   c = 0.5_dp
   m1 = Edweibull3(c,0.0_dp)
   m2 = E2dweibull3(c,0.0_dp)
   if (abs(m1-exp(-c)/(1.0_dp-exp(-c))) > 1.0e-14_dp) fails=fails+1
   if (abs(m2-exp(-c)*(1.0_dp+exp(-c))/(1.0_dp-exp(-c))**2) > 1.0e-14_dp) &
      fails=fails+1

   ! Independent direct-sum check for a nontrivial type-I case.
   q = 0.55_dp
   m1 = Edweibull(q,1.4_dp,eps=1.0e-8_dp,nmax=100000)
   if (abs(m1-direct_mean(q,1.4_dp)) > 2.0e-8_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_moments: FAIL", fails
      error stop 1
   end if
   print *, "test_moments: PASS"

contains
   real(dp) function direct_mean(q,beta) result(m)
      real(dp), intent(in) :: q,beta
      integer(i64) :: k
      real(dp) :: term
      m = 0.0_dp
      do k = 1_i64, 100000_i64
         term = ddweibull(k,q,beta,.false.)
         m = m+real(k,dp)*term
         if (term < 1.0e-16_dp) exit
      end do
   end function direct_mean
end program test_moments
