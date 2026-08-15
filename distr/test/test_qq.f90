program test_qq
   use distr_kinds, only : dp
   use distr_core, only : distribution_t, normal_dist, chisq_dist, poisson_dist
   use distr_qq, only : qqbounds, qq_ks_critical
   implicit none
   integer :: fails
   type(distribution_t) :: d
   real(dp) :: x(1),b(1,4),v

   fails=0

   ! Upstream distr Examples fixture (alpha=.95, n=30), with x recovered as
   ! the midpoint of the printed symmetric pointwise interval.
   d=normal_dist()
   x(1)=-1.91318805_dp
   call qqbounds(x,d,0.95_dp,30,b,exact_pointwise=.true.,exact_simultaneous=.true.)
   call check_close('exact KS critical',qq_ks_critical(0.95_dp,30,.true.),1.32386443_dp,2.0e-8_dp)
   call check_close('exact normal sim left',b(1,1),-4.7534243_dp,5.0e-7_dp)
   call check_close('exact normal sim right',b(1,2),-0.6140964_dp,5.0e-5_dp)
   call check_close('exact normal pw left',b(1,3),-2.9401996_dp,5.0e-7_dp)
   call check_close('exact normal pw right',b(1,4),-0.8861765_dp,5.0e-7_dp)

   ! Upstream asymptotic qqbounds behavior is intentionally reproduced,
   ! including the complement convention in .q2kolmogorov.
   d=chisq_dist(4.0_dp)
   x(1)=0.6476086_dp
   call qqbounds(x,d,0.95_dp,30,b,exact_pointwise=.false.,exact_simultaneous=.false.)
   call check_close('asym KS critical',qq_ks_critical(0.95_dp,30,.false.),0.51961038_dp,2.0e-8_dp)
   call check_close('asym chisq sim left',b(1,1),0.00282976_dp,2.0e-8_dp)
   call check_close('asym chisq sim right',b(1,2),1.29183_dp,2.0e-5_dp)
   call check_close('asym chisq pw left',b(1,3),0.0321572_dp,2.0e-6_dp)
   call check_close('asym chisq pw right',b(1,4),1.26306_dp,2.0e-6_dp)

   d=poisson_dist(3.0_dp)
   v=d%cdf_left(2.0_dp)
   call check_close('left CDF at atom',v,d%cdf(1.0_dp),2.0e-14_dp)

   if (fails==0) then
      print '(a)', 'test_qq: PASS'
   else
      print '(a,i0)', 'test_qq: FAIL ',fails
      error stop 1
   end if

contains

   subroutine check_close(name,got,want,tol)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: got,want,tol
      if (abs(got-want)>tol) then
         fails=fails+1
         print '(a,2es24.15)', trim(name)//' got/want: ',got,want
      end if
   end subroutine check_close

end program test_qq
