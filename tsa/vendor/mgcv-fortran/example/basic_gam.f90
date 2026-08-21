program basic_gam
   use mgcv
   implicit none
   integer, parameter :: n = 160, k = 12
   real(dp) :: x(n), y(n)
   real(dp), allocatable :: basis(:, :), design(:, :), penalties(:, :, :), full_penalty(:, :)
   type(smooth_spec_t) :: smooth
   type(gam_model_t) :: model
   type(family_t) :: family
   integer :: i, status

   do i = 1, n
      x(i) = real(i - 1, dp) / real(n - 1, dp)
      y(i) = sin(2.0_dp * pi_dp * x(i)) + 0.35_dp * cos(6.0_dp * pi_dp * x(i)) + &
             0.08_dp * sin(37.0_dp * x(i))
   end do

   call construct_ps_smooth(x, k, basis, smooth, status)
   if (status /= 0) error stop 'unable to construct P-spline basis'

   design = append_columns(reshape([(1.0_dp, i=1,n)], [n,1]), basis)
   allocate(penalties(k + 1, k + 1, 1)); penalties = 0.0_dp
   full_penalty = embed_penalty(smooth%penalties(:, :, 1), 2, k + 1)
   penalties(:, :, 1) = full_penalty

   family%id = family_gaussian
   call gam_fit(design, y, penalties, model, status, family=family, method=method_gcv)
   if (status /= 0) error stop 'GAM fit failed'

   call model%summary()
   write(*,'(a,f10.6)') 'RMSE against observations: ', &
      sqrt(sum((model%fitted - y)**2) / real(n, dp))
end program basic_gam
