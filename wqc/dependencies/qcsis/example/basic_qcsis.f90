program basic_qcsis
   use qcsis_mod, only : cqcsis, dp, qcsis, screening_result
   implicit none

   integer, parameter :: n = 100
   integer, parameter :: p = 12
   real(dp) :: x(n, p), y(n), noise(n)
   type(screening_result) :: fit
   integer, allocatable :: seed(:)
   integer :: i, seed_size, stat
   character(len=:), allocatable :: errmsg

   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   seed = [(104729 + 37 * i, i = 1, seed_size)]
   call random_seed(put=seed)

   call random_number(x)
   call random_number(noise)
   x = 2.0_dp * x - 1.0_dp
   noise = 0.1_dp * (2.0_dp * noise - 1.0_dp)
   y = 3.0_dp * x(:, 1) - 2.0_dp * x(:, 2) + 1.5_dp * x(:, 3) + noise

   fit = qcsis(x, y, 3, stat, errmsg)
   if (stat /= 0) error stop errmsg

   print '(a)', "QCSIS selected predictor indices:"
   print '(*(i0,1x))', fit%selected

   fit = cqcsis(x, y, 3, stat, errmsg)
   if (stat /= 0) error stop errmsg

   print '(a)', "CQCSIS selected predictor indices:"
   print '(*(i0,1x))', fit%selected
   print '(a)', "CQCSIS weights:"
   do i = 1, p
      print '(i3,2x,es14.6)', i, fit%w(i)
   end do
end program basic_qcsis
