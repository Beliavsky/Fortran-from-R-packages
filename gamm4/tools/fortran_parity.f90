program fortran_parity
   use gamm4
   implicit none
   integer, parameter :: n = 18
   type(gamm4_smooth_t) :: smooths(1)
   type(random_term_t) :: terms(1)
   type(gamm4_result_t) :: fit
   type(gamm4_control_t) :: control
   real(dp) :: x(n), y(n), fixed(n, 1), ge, err
   integer :: i, status

   do i = 1, n
      x(i) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
      ge = merge(-0.28_dp, 0.28_dp, i <= n / 2)
      err = 0.025_dp * sin(1.31_dp * real(i, dp))
      y(i) = 0.9_dp + 1.1_dp * x(i) + 0.75_dp * x(i) * x(i) - 0.55_dp * x(i) ** 3 + ge + err
   end do
   fixed(:, 1) = 1.0_dp
   allocate(smooths(1)%basis(n, 3), smooths(1)%spec%penalties(3, 3, 1))
   smooths(1)%basis(:, 1) = x
   smooths(1)%basis(:, 2) = x * x
   smooths(1)%basis(:, 3) = x ** 3
   smooths(1)%spec%penalties = 0.0_dp
   smooths(1)%spec%penalties(2, 2, 1) = 1.0_dp
   smooths(1)%spec%penalties(3, 3, 1) = 1.0_dp
   smooths(1)%label = 's(x)'
   allocate(terms(1)%z(n, 1), terms(1)%group(n))
   terms(1)%z = 1.0_dp
   terms(1)%group(1:n / 2) = 1
   terms(1)%group(n / 2 + 1:n) = 2
   terms(1)%n_levels = 2
   terms(1)%covariance_structure = covariance_diagonal
   terms(1)%name = 'group'
   control%max_outer = 8
   control%tolerance = 1.0e-5_dp
   call gamm4_fit(y, fixed, smooths, fit, status, random_terms=terms, reml=.false., control=control)
   if (status /= 0) error stop 'parity fit failed'
   open(unit=10, file='fortran_parity.csv', status='replace', action='write')
   write(10, '(a)') 'name,value'
   write(10, '(a,es24.16)') 'objective,', fit%deviance
   write(10, '(a,es24.16)') 'scale,', fit%scale
   write(10, '(a,es24.16)') 'sp,', fit%sp(1)
   do i = 1, size(fit%coefficients)
      write(10, '(a,i0,a,es24.16)') 'coef', i, ',', fit%coefficients(i)
   end do
   do i = 1, size(fit%variance_parameters)
      write(10, '(a,i0,a,es24.16)') 'eta', i, ',', fit%variance_parameters(i)
   end do
   do i = 1, size(fit%covariance, 1)
      write(10, '(a,i0,a,es24.16)') 'vcov_diag', i, ',', fit%covariance(i, i)
   end do
   close(10)
end program fortran_parity
