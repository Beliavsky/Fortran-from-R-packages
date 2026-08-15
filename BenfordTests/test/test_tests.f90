program test_tests
   use benfordtests
   implicit none
   real(dp), allocatable :: x(:),h0(:)
   type(benford_test_result_t) :: r
   type(jointdigit_result_t) :: jr
   type(significant_digit_analysis_t) :: ar
   integer :: fails,i,j,k
   real(dp),parameter::tol=2.0d-12
   fails=0
   allocate(x(12));k=0
   do i=1,9
      do j=1,merge(3,merge(2,1,i==2),i==1)
         k=k+1;x(k)=real(i,dp)
      end do
   end do
   call chisq_benftest(x,1,r)
   call check(abs(r%statistic-1.0958006038990171_dp)<tol,'chi-square statistic')
   call check(abs(r%p_value-0.9975676371883492_dp)<2.0d-12,'chi-square p-value')
   call ks_benftest(x,1,r,pvalsims=100)
   call check(abs(r%statistic-0.4005771787895563_dp)<tol,'KS statistic')
   call mdist_benftest(x,1,r,pvalsims=100)
   call check(abs(r%statistic-0.17677309040006986_dp)<tol,'m-distance statistic')
   call edist_benftest(x,1,r,pvalsims=100)
   call check(abs(r%statistic-0.30952309639746184_dp)<tol,'e-distance statistic')
   call usq_benftest(x,1,r,pvalsims=100)
   call check(abs(r%statistic-0.01579298003346328_dp)<tol,'U-square statistic')
   call meandigit_benftest(x,1,r)
   call check(abs(r%statistic-0.11566974390945735_dp)<tol,'mean-digit statistic')
   call jpsq_benftest(x,1,r,pvalsims=100)
   call check(abs(r%statistic-0.9102630309457205_dp)<tol,'JP-square statistic')
   call jointdigit_benftest(x,1,jr,'all')
   call check(abs(jr%statistic-1.0958006038990171_dp)<2.0d-10,'joint all statistic')
   call check(jr%df==8,'joint all df')
   call jointdigit_benftest(x,1,jr,'kaiser')
   call check(jr%df>=1 .and. jr%df<=8,'joint kaiser df')
   call jointdigit_benftest_indices(x,1,[1,2],jr)
   call check(jr%df==2 .and. all(jr%eigenvalues_tested==[1,2]),'joint explicit eigen indices')
   call simulate_h0('chisq',50,1,200,h0)
   call check(size(h0)==200 .and. all(h0>=0.0_dp),'simulate H0')
   call signifd_analysis(x,1,ar,alpha=[0.05_dp],relative=.true.)
   call check(size(ar%frequency)==9 .and. size(ar%confidence,1)==3,'digit analysis shapes')
   call check(abs(sum(ar%frequency)-1.0_dp)<tol,'digit analysis frequencies')
   if(fails==0)then
      print '(a)','test_tests: PASS'
   else
      print '(a,i0)','test_tests: FAIL ',fails
      error stop 1
   end if
contains
   subroutine check(ok,name)
      logical,intent(in)::ok;character(len=*),intent(in)::name
      if(.not.ok)then;fails=fails+1;print '(a,a)','FAIL: ',trim(name);end if
   end subroutine
end program
