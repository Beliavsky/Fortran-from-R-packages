program compare_rfortran_core
   use r_kinds, only : dp
   use r_rolling, only : r_roll_mean_right, r_roll_mean_valid
   use r_special, only : r_digamma, r_log_beta, r_trigamma
   implicit none

   real(dp), parameter :: points(6) = [0.1_dp, 0.5_dp, 1.0_dp, 2.5_dp, 10.0_dp, 100.0_dp]
   real(dp) :: x(100), value, t0, t1
   real(dp), allocatable :: rolling(:)
   integer :: i, j, reps, unit
   character(512) :: output

   call get_command_argument(1, output)
   if (len_trim(output) == 0) output = 'fortran_results.csv'
   x = [(real(i, dp), i=1,100)]
   reps = 20000

   open(newunit=unit, file=trim(output), status='replace')
   write(unit, '(a)') 'case,value,seconds,abs_tol,rel_tol'

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_digamma(points)*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('digamma_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_trigamma(points)*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('trigamma_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      value = sum(r_log_beta(points, points(size(points):1:-1))*[(real(j, dp), j=1,size(points))])
   end do
   call cpu_time(t1)
   call emit('log_beta_checksum', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_mean_valid(x, 5, rolling)
      value = sum(rolling*[(real(j, dp), j=1,size(rolling))])
   end do
   call cpu_time(t1)
   call emit('rolling_mean_valid', value, t1 - t0)

   call cpu_time(t0)
   do i = 1, reps
      call r_roll_mean_right(x, 5, rolling)
      value = sum(rolling(5:)*[(real(j, dp), j=5,100)])
   end do
   call cpu_time(t1)
   call emit('rolling_mean_right', value, t1 - t0)

   close(unit)

contains

   subroutine emit(name, result, seconds)
      character(*), intent(in) :: name
      real(dp), intent(in) :: result, seconds

      write(unit, '(a,",",es25.16e3,",",es16.8,",",es12.4,",",es12.4)') &
         trim(name), result, seconds, 1.0e-11_dp, 1.0e-11_dp
   end subroutine emit

end program compare_rfortran_core
