program mbend_example
   use mbend_api, only : bend, bend_result, dp, mbend_success
   implicit none

   real(dp) :: v(5, 5)
   type(bend_result) :: fit
   integer :: info

   v = reshape([ &
      100.0_dp, 95.0_dp, 80.0_dp, 40.0_dp, 40.0_dp, &
       95.0_dp,100.0_dp, 95.0_dp, 80.0_dp, 40.0_dp, &
       80.0_dp, 95.0_dp,100.0_dp, 95.0_dp, 80.0_dp, &
       40.0_dp, 80.0_dp, 95.0_dp,100.0_dp, 95.0_dp, &
       40.0_dp, 40.0_dp, 80.0_dp, 95.0_dp,100.0_dp], [5, 5])

   call bend(v, fit, info)
   if (info /= mbend_success) error stop 'mbend example failed'

   print '(a,i0)', 'iterations = ', fit%iterations
   print '(a,es14.6)', 'minimum final eigenvalue = ', minval(fit%final_ev)
   print '(a)', 'bent matrix:'
   call print_matrix(fit%bent)

contains

   subroutine print_matrix(a)
      real(dp), intent(in) :: a(:, :) !! Matrix printed one row per output record.

      integer :: i

      do i = 1, size(a, 1)
         print '(*(f12.6,1x))', a(i, :)
      end do
   end subroutine print_matrix

end program mbend_example
