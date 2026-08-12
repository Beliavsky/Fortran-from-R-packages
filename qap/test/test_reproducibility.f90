program test_reproducibility
   use qap
   implicit none
   integer, parameter :: n = 9
   real(dp) :: A(n,n), B(n,n)
   integer :: i, j
   type(qap_control_t) :: ctl
   type(qap_result_t) :: r1, r2

   do j = 1, n
      do i = 1, n
         A(i,j) = real((i-j)*(i-j),dp)
         B(i,j) = real(mod(i*j + 3*i + 2*j,17),dp)
      end do
   end do
   A = 0.5_dp * (A + transpose(A))
   B = 0.5_dp * (B + transpose(B))

   ctl%rep = 4
   ctl%seed = 987654_i64
   call qap_solve(A,B,r1,ctl)
   call qap_solve(A,B,r2,ctl)
   if (any(r1%permutation /= r2%permutation)) error stop 'seed is not reproducible'
   if (abs(r1%objective-r2%objective) > 0.0_dp) error stop 'objective is not reproducible'
   if (r1%attempted_swaps /= r2%attempted_swaps) error stop 'trajectory differs'
   print *, 'test_reproducibility: PASS'
end program test_reproducibility
