program basic_kde1d
   use kde1d_api, only : dkde1d, dp, kde1d_fit, kde1d_model, pkde1d, qkde1d
   implicit none

   type(kde1d_model) :: model
   real(dp) :: x(40)
   real(dp) :: points(5)
   integer :: i
   integer :: ierr

   do i = 1, size(x)
      x(i) = sin(0.22_dp * real(i, dp)) + 0.05_dp * real(i, dp)
   end do
   call kde1d_fit(x, model, ierr, deg=2, boundary_repair=.false.)
   if (ierr /= 0) error stop 'kde1d_fit failed'

   points = [-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, 2.0_dp]
   print '(a,5f10.4)', 'x:   ', points
   print '(a,5f10.4)', 'pdf: ', dkde1d(points, model)
   print '(a,5f10.4)', 'cdf: ', pkde1d(points, model)
   print '(a,5f10.4)', 'q:   ', qkde1d([0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.9_dp], model)
   print '(a,f10.5)', 'bandwidth: ', model%bw
end program basic_kde1d
