program test_functions
   use fgarch_kinds, only : dp
   use tvgarchkf
   use test_support
   implicit none
   type(tv_function_spec) :: spec
   real(dp), allocatable :: values(:)
   real(dp) :: expected(4), u
   integer :: i, status
   character(len=160) :: message

   spec = make_tv_function([1.0_dp,2.0_dp,3.0_dp],tv_polynomial)
   values = evaluate_tv_function(spec,4,status,message)
   do i = 1, 4
      u = real(i,dp)/4.0_dp
      expected(i) = 1.0_dp+2.0_dp*u+3.0_dp*u*u
   end do
   call assert_true(status == 0,'polynomial status')
   call assert_all_close(values,expected,1.0e-13_dp,'polynomial values')

   spec = make_tv_function([1.0_dp,2.0_dp],tv_nonlinear,[0.0_dp,0.5_dp])
   values = evaluate_tv_function(spec,4,status,message)
   do i = 1, 4
      u = real(i,dp)/4.0_dp
      expected(i) = 1.0_dp+2.0_dp*sqrt(u)
   end do
   call assert_all_close(values,expected,1.0e-13_dp,'nonlinear values')

   spec = make_tv_function([1.0_dp,2.0_dp,3.0_dp],tv_trigonometric, &
                           trig_kind=trig_cosine,argument_kind=arg_identity)
   values = evaluate_tv_function(spec,4,status,message)
   do i = 1, 4
      u = real(i,dp)/4.0_dp
      expected(i) = 1.0_dp+5.0_dp*cos(u)
   end do
   call assert_all_close(values,expected,1.0e-13_dp,'trigonometric values')

   spec = make_tv_function([0.5_dp,0.25_dp],tv_trigonometric, &
                           trig_kind=trig_sine,argument_kind=arg_three_one_minus_log)
   values = evaluate_tv_function(spec,4,status,message)
   do i = 1, 4
      u = real(i,dp)/4.0_dp
      expected(i) = 0.5_dp+0.25_dp*sin(3.0_dp*(1.0_dp-log(u)))
   end do
   call assert_all_close(values,expected,1.0e-13_dp,'documented argument values')

   spec = make_tv_function([1.0_dp],99)
   values = evaluate_tv_function(spec,4,status,message)
   call assert_true(status /= 0,'invalid function kind rejected')
   write(*,'(a)') 'test_functions: PASS'
end program test_functions
