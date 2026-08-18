program test_friedman_tables
   use suppdists
   implicit none
   integer :: fails
   fails=0
   call check(pfriedman(1.0_dp,3,2),0.5_dp,2e-12_dp,'pfriedman exact')
   call check(dfriedman(1.0_dp,3,2),1.0_dp/3.0_dp,2e-12_dp,'dfriedman exact')
   call check(pspearman(-0.5_dp,3),0.5_dp,2e-12_dp,'pspearman exact')
   call check(qspearman(0.5_dp,3),-0.5_dp,2e-12_dp,'qspearman exact')
   if(fails==0)then
      print '(a)','test_friedman_tables: PASS'
   else
      error stop 1
   end if
contains
   subroutine check(got,want,tol,name)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::name
      if(abs(got-want)>tol)then
         print '(a,2es24.14)',trim(name)//' FAIL ',got,want
         fails=fails+1
      end if
   end subroutine check
end program test_friedman_tables
