! SPDX-License-Identifier: GPL-2.0-or-later
program fit_mvnormal_mixture
  use, intrinsic :: iso_fortran_env, only : int64
  use mixtools, only : dp, em_control, mv_mixture_result, mvnormalmixEM, &
    MIXTOOLS_SUCCESS, mixtools_status_message
  implicit none

  integer, parameter :: component_count = 2
  character(len=512) :: data_file, report_file
  integer :: data_unit, report_unit, ios, n, p, i
  integer(int64) :: clock_rate, overall_start, read_start, read_end
  integer(int64) :: fit_start, fit_end
  real(dp) :: read_seconds, fit_seconds, overall_seconds
  real(dp), allocatable :: observations(:,:)
  type(em_control) :: control
  type(mv_mixture_result) :: fit

  call system_clock(count=overall_start, count_rate=clock_rate)
  if (clock_rate <= 0_int64) error stop 'System clock is unavailable.'
  data_file = 'example/mvnormal_mixture_data.txt'
  report_file = 'example/mvnormal_mixture_fit_fortran.txt'
  if (command_argument_count() >= 1) call get_command_argument(1, data_file)
  if (command_argument_count() >= 2) call get_command_argument(2, report_file)

  call system_clock(count=read_start)
  open(newunit=data_unit, file=trim(data_file), status='old', action='read', &
    iostat=ios)
  if (ios /= 0) error stop 'Could not open the observation file.'

  read(data_unit, *, iostat=ios) n, p
  if (ios /= 0 .or. n < component_count .or. p < 1) &
    error stop 'Invalid observation-file header.'

  allocate(observations(n,p))
  do i = 1, n
    read(data_unit, *, iostat=ios) observations(i,:)
    if (ios /= 0) error stop 'Invalid or incomplete observation data.'
  end do
  close(data_unit)
  call system_clock(count=read_end)
  read_seconds = elapsed_seconds(read_start, read_end, clock_rate)

  control%tolerance = 1.0e-8_dp
  control%max_iterations = 1000
  control%ridge = 1.0e-8_dp
  call system_clock(count=fit_start)
  call mvnormalmixEM(observations, component_count, fit, control)
  call system_clock(count=fit_end)
  fit_seconds = elapsed_seconds(fit_start, fit_end, clock_rate)
  overall_seconds = elapsed_seconds(overall_start, fit_end, clock_rate)

  if (fit%status /= MIXTOOLS_SUCCESS) then
    write(*, '(a)') 'Fit failed: ' // mixtools_status_message(fit%status)
    error stop 1
  end if

  call write_report(6, fit, n, p, read_seconds, fit_seconds, overall_seconds)
  open(newunit=report_unit, file=trim(report_file), status='replace', &
    action='write', iostat=ios)
  if (ios /= 0) error stop 'Could not open the Fortran report file.'
  call write_report(report_unit, fit, n, p, read_seconds, fit_seconds, &
    overall_seconds)
  close(report_unit)
  write(*, '(a)') 'Wrote fit report to ' // trim(report_file)

contains

  pure function elapsed_seconds(start_count, end_count, rate) result(seconds)
    integer(int64), intent(in) :: start_count, end_count, rate
    real(dp) :: seconds

    seconds = real(end_count - start_count, dp) / real(rate, dp)
  end function elapsed_seconds

  subroutine write_report(unit, result, observation_count, dimension_count, &
      input_seconds, mixture_seconds, total_seconds)
    integer, intent(in) :: unit, observation_count, dimension_count
    type(mv_mixture_result), intent(in) :: result
    real(dp), intent(in) :: input_seconds, mixture_seconds, total_seconds
    integer :: component, position, row
    integer :: order(component_count)

    order = [1, 2]
    if (result%mu(1,1) > result%mu(1,2)) order = [2, 1]

    write(unit, '(a)') 'Fortran mixtools mvnormalmixEM fit'
    write(unit, '(a,i0,a,i0)') 'observations: ', observation_count, &
      ' dimensions: ', dimension_count
    write(unit, '(a,es20.12)') 'loglik: ', result%loglik
    write(unit, '(a,i0)') 'iterations: ', result%iterations
    write(unit, '(a,l1)') 'converged: ', result%converged
    write(unit, '(a,f0.6)') 'read seconds: ', input_seconds
    write(unit, '(a,f0.6)') 'fit seconds: ', mixture_seconds
    write(unit, '(a,f0.6)') 'overall seconds: ', total_seconds
    do position = 1, component_count
      component = order(position)
      write(unit, '(a,i0)') 'component ', position
      write(unit, '(a,es20.12)') '  weight: ', result%lambda(component)
      write(unit, '(a,*(es20.12,1x))') '  mean: ', result%mu(:,component)
      write(unit, '(a)') '  covariance:'
      do row = 1, dimension_count
        write(unit, '(*(es20.12,1x))') result%sigma(row,:,component)
      end do
    end do
  end subroutine write_report

end program fit_mvnormal_mixture
