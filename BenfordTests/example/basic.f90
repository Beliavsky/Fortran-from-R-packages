program basic
   use benfordtests
   implicit none
   real(dp), allocatable :: x(:)
   type(benford_test_result_t) :: result
   allocate(x(5000));call rbenf(x)
   call chisq_benftest(x,1,result)
   print '(a,f10.5)','Chi-square statistic: ',result%statistic
   print '(a,f10.5)','Asymptotic p-value:   ',result%p_value
   print '(a)','First-digit Benford probabilities:'
   print '(9f8.5)',pbenf(1)
end program
