program geepack_example
   use geepack
   implicit none

   real(dp) :: y(18)
   real(dp) :: x(18, 2)
   integer :: cluster_sizes(6)
   real(dp) :: beta_initial(2)
   type(gee_spec) :: spec
   type(gee_result) :: fit
   integer :: cluster
   integer :: visit
   integer :: row

   cluster_sizes = 3
   x(:, 1) = 1.0_dp
   row = 0
   do cluster = 1, 6
      do visit = 1, 3
         row = row + 1
         x(row, 2) = real(visit - 2, dp)
         y(row) = 2.0_dp + 0.7_dp * x(row, 2) + 0.12_dp * real(cluster - 3, dp)
      end do
   end do

   beta_initial = [2.0_dp, 0.5_dp]
   spec%corstr = COR_EXCHANGEABLE
   spec%scale_fixed = .true.
   spec%scale_value = 1.0_dp
   spec%tolerance = 1.0e-10_dp
   allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
   spec%mean_links = LINK_IDENTITY
   spec%variance_codes = VAR_GAUSSIAN
   spec%scale_links = LINK_IDENTITY

   call fit_geese(y, x, cluster_sizes, beta_initial, spec, fit)
   if (fit%error /= GEE_OK) error stop "geepack example fit failed"

   print '(a,2f12.6)', 'beta:  ', fit%beta
   print '(a,f12.6)', 'alpha: ', fit%alpha(1)
end program geepack_example
