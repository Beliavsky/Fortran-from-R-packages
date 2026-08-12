program test_swap_delta
   use qap
   implicit none
   integer, parameter :: n = 6
   real(dp) :: A(n,n), B(n,n), z0, z1, delta
   integer :: perm(n), trial(n), i, j, k, tmp

   do j = 1, n
      do i = 1, n
         A(i,j) = real(abs(i-j) + mod(i+j,3), dp)
         B(i,j) = real(abs(i-j) + mod(i*j,5), dp)
      end do
   end do
   A = 0.5_dp * (A + transpose(A))
   B = 0.5_dp * (B + transpose(B))
   perm = [3, 6, 1, 5, 2, 4]
   z0 = qap_obj(A, B, perm)

   do i = 1, n-1
      do j = i+1, n
         trial = perm
         tmp = trial(i)
         trial(i) = trial(j)
         trial(j) = tmp
         delta = qap_swap_delta(A, B, perm, i, j)
         z1 = qap_obj(A, B, trial)
         if (abs((z1-z0) - delta) > 1.0e-10_dp) then
            write(*,*) i, j, z1-z0, delta
            error stop 'swap delta mismatch'
         end if
      end do
   end do

   k = 1
   if (k /= 1) error stop 'unreachable'
   print *, 'test_swap_delta: PASS'
end program test_swap_delta
