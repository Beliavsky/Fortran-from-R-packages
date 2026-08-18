program test_semidefinite
   use misc_tools
   implicit none
   real(dp) :: a(3,3),psd(3,3),sing(4,4),neg(4,4),zero(3,3),indef(2,2)
   real(dp), allocatable :: eig(:)
   integer :: fails,status

   fails = 0

   a = reshape([1.0_dp,2.0_dp,3.0_dp, &
                4.0_dp,5.0_dp,6.0_dp, &
                7.0_dp,8.0_dp,10.0_dp],[3,3])
   psd = matmul(transpose(a),a)

   if (.not. semidefiniteness(psd,method="det")) fails=fails+1
   if (.not. semidefiniteness(psd,method="eigen")) fails=fails+1
   if (semidefiniteness(psd,positive=.false.,method="eigen")) fails=fails+1

   sing = 0.0_dp
   sing(1:3,1:3) = psd
   sing(1:3,4) = -sum(psd,dim=2)
   sing(4,1:3) = -sum(psd,dim=1)
   sing(4,4) = sum(psd)
   if (.not. semidefiniteness(sing,method="eigen")) fails=fails+1

   neg = 0.0_dp
   neg(1,1)=-1.0_dp; neg(2,2)=-1.0_dp; neg(3,3)=-1.0_dp; neg(4,4)=-1.0_dp
   if (semidefiniteness(neg,method="eigen")) fails=fails+1
   if (.not. semidefiniteness(neg,positive=.false.,method="eigen")) fails=fails+1

   zero = 0.0_dp
   if (.not. semidefiniteness(zero,method="det")) fails=fails+1
   if (.not. semidefiniteness(zero,positive=.false.,method="eigen")) fails=fails+1

   indef = reshape([1.0_dp,0.0_dp,0.0_dp,-1.0_dp],[2,2])
   if (semidefiniteness(indef,method="det")) fails=fails+1
   if (semidefiniteness(indef,positive=.false.,method="eigen")) fails=fails+1

   call symmetric_eigenvalues(psd,eig,status)
   if (status /= 0 .or. minval(eig) < -1.0e-10_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_semidefinite: FAIL", fails
      error stop 1
   end if
   print *, "test_semidefinite: PASS"
end program test_semidefinite
