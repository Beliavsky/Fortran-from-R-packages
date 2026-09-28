program pdfcluster_example
   use pdfcluster, only : dp, groups_at_stage, pdfcluster_fit, pdfcluster_result
   implicit none

   type(pdfcluster_result) :: fit
   real(dp) :: x(10, 1)
   integer, allocatable :: groups(:)
   logical :: ok

   x(:, 1) = [-5.0_dp, -4.9_dp, -4.8_dp, -4.7_dp, -4.6_dp, &
                4.6_dp, 4.7_dp, 4.8_dp, 4.9_dp, 5.0_dp]
   call pdfcluster_fit(x, fit, ok, n_grid=10, n_stage=5)
   if (.not. ok) error stop "pdfCluster example fit failed"
   call groups_at_stage(fit, groups=groups, success=ok)
   if (.not. ok) error stop "pdfCluster example group extraction failed"

   print *, "Detected cluster cores:", fit%noc
   print *, "Final groups:", groups
end program pdfcluster_example
