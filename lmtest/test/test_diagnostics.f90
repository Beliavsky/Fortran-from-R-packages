program test_diagnostics
   use lmtest, only : dp, test_result, bg_test_result, &
      breusch_godfrey_test, breusch_pagan_test, goldfeld_quandt_test, &
      harvey_collier_test, rainbow_test, reset_test, durbin_watson_test, &
      harrison_mccabe_test
   implicit none
   integer, parameter :: n = 80
   real(dp) :: x(n,3), y(n)
   type(bg_test_result) :: bg
   type(test_result) :: tr

   call make_data(x,y)

   bg = breusch_godfrey_test(x,y,order=2)
   call check(bg%test%statistic, 75.54184783704888_dp, 3.0e-10_dp)
   call check(bg%test%p_value, 3.947263927716519e-17_dp, 2.0e-28_dp)

   bg = breusch_godfrey_test(x,y,order=2,use_f=.true.)
   call check(bg%test%statistic, 635.4245414571315_dp, 3.0e-9_dp)

   tr = breusch_pagan_test(x,y,x,studentize=.true.)
   call check(tr%statistic, 7.271207977024421_dp, 3.0e-11_dp)
   call check(tr%p_value, 0.02636800360563612_dp, 3.0e-13_dp)

   tr = breusch_pagan_test(x,y,x,studentize=.false.)
   call check(tr%statistic, 3.0794322991952825_dp, 3.0e-11_dp)

   tr = goldfeld_quandt_test(x,y)
   call check(tr%statistic, 1.641760171679299_dp, 2.0e-12_dp)
   call check(tr%p_value, 0.06809147841060445_dp, 3.0e-13_dp)

   tr = harvey_collier_test(x,y)
   call check(tr%statistic, 1.4066348820050234_dp, 2.0e-8_dp)
   call check(tr%p_value, 0.16361172096894164_dp, 2.0e-9_dp)

   tr = rainbow_test(x,y)
   call check(tr%statistic, 1.130726261332067_dp, 2.0e-12_dp)
   call check(tr%p_value, 0.3542824750045829_dp, 3.0e-13_dp)

   tr = reset_test(x,y,powers=[2])
   call check(tr%statistic, 1.7331269669918195_dp, 3.0e-10_dp)
   call check(tr%p_value, 0.19196848252011878_dp, 3.0e-12_dp)

   tr = durbin_watson_test(x,y,exact=.false.)
   call check(tr%statistic, 0.7049015438539227_dp, 2.0e-13_dp)
   if (tr%p_value <= 0.0_dp .or. tr%p_value >= 1.0_dp) error stop 1

   tr = durbin_watson_test(x,y,exact=.true.)
   call check(tr%statistic, 0.7049015438539227_dp, 2.0e-13_dp)
   if (tr%p_value <= 0.0_dp .or. tr%p_value >= 1.0_dp) error stop 1

   tr = harrison_mccabe_test(x,y,nsim=0)
   if (tr%statistic <= 0.0_dp .or. tr%statistic >= 1.0_dp) error stop 1

contains
   subroutine make_data(a,b)
      real(dp), intent(out) :: a(:,:), b(:)
      integer :: j
      real(dp) :: v
      do j=1,size(b)
         v=real(j,dp)
         a(j,1)=1.0_dp
         a(j,2)=(v-40.0_dp)/20.0_dp
         a(j,3)=sin(0.3_dp*v)
         b(j)=1.0_dp+2.0_dp*a(j,2)-0.7_dp*a(j,3)+ &
            (0.18_dp+0.002_dp*v)*sin(0.91_dp*v)+0.08_dp*cos(0.17_dp*v)
      end do
   end subroutine make_data
   subroutine check(actual, expected, tol)
      real(dp), intent(in) :: actual, expected, tol
      if (abs(actual-expected)>tol) then
         print *, actual, expected, tol
         error stop 1
      end if
   end subroutine check
end program test_diagnostics
