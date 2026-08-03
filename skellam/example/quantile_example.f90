program quantile_example
   use skellam, only : dp, i8, qskellam
   implicit none
   real(dp), parameter :: p(5) = [0.01_dp, 0.05_dp, 0.5_dp, 0.95_dp, 0.99_dp]
   integer(i8) :: q(5)
   integer :: i

   q = qskellam(p, 3.0_dp, 4.0_dp)
   do i = 1, size(p)
      print '(f8.3,2x,i0)', p(i), q(i)
   end do
end program quantile_example
