program test_model
   use highs
   implicit none
   type(highs_model) :: model
   integer(highs_int) :: status
   real(dp) :: a(2,2), q(2,2)

   a = reshape([1.0_dp,3.0_dp,2.0_dp,4.0_dp], [2,2])
   q = reshape([2.0_dp,0.5_dp,0.5_dp,4.0_dp], [2,2])
   call highs_model_from_dense(model, [1.0_dp,2.0_dp], [0.0_dp,0.0_dp], &
      [10.0_dp,10.0_dp], status, a=a, lhs=[1.0_dp,2.0_dp], &
      rhs=[3.0_dp,4.0_dp], vartype=[highs_var_continuous,highs_var_integer], &
      maximum=.true., offset=7.0_dp, q=q)
   if (status /= highs_status_ok) error stop "model construction failed"
   if (.not. model%valid()) error stop "model validation failed"
   if (model%sense /= highs_maximize .or. .not. model%has_hessian) error stop "model fields wrong"
   if (model%a%nnz() /= 4 .or. model%q%nnz() /= 3) error stop "sparse model counts wrong"
   print *, "test_model: PASS"
end program test_model
