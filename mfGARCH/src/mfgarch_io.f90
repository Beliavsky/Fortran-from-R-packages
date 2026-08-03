! SPDX-License-Identifier: MIT
module mfgarch_io
  use mfgarch_kinds, only : dp
  use mfgarch_components, only : model_parameters, parameter_names
  use mfgarch_math, only : finite_value
  use mfgarch_status, only : mfgarch_success
  use mfgarch_types, only : mfgarch_fit_result, mfgarch_simulation
  implicit none
  private

  public :: print_fit_summary, write_simulation_csv

contains

  subroutine print_fit_summary(result, unit)
    type(mfgarch_fit_result), intent(in) :: result
    integer, intent(in), optional :: unit
    real(dp), allocatable :: parameters(:)
    character(len=16), allocatable :: names(:)
    integer :: output_unit, i

    output_unit = 6
    if (present(unit)) output_unit = unit
    call model_parameters(result%model, parameters)
    call parameter_names(result%model, names)
    write(output_unit,'(a)') 'mfGARCH fit'
    write(output_unit,'(a,es16.8)') 'log likelihood: ', result%log_likelihood
    write(output_unit,'(a,es16.8)') 'BIC:            ', result%bic
    write(output_unit,'(a,i0)') 'status:         ', result%status
    write(output_unit,'(a,a)') 'message:        ', trim(result%message)
    write(output_unit,'(a)') ''
    write(output_unit,'(a16,2x,a16,2x,a16)') 'parameter', 'estimate', 'robust s.e.'
    do i = 1, size(parameters)
      if (allocated(result%robust_standard_error) .and. &
          size(result%robust_standard_error) == size(parameters) .and. &
          finite_value(result%robust_standard_error(i))) then
        write(output_unit,'(a16,2x,es16.8,2x,es16.8)') trim(names(i)), parameters(i), &
          result%robust_standard_error(i)
      else
        write(output_unit,'(a16,2x,es16.8,2x,a16)') trim(names(i)), parameters(i), 'NA'
      end if
    end do
    if (result%model%k > 0) then
      write(output_unit,'(a,es16.8)') 'variance ratio: ', result%variance_ratio
      write(output_unit,'(a,es16.8)') 'tau forecast:   ', result%tau_forecast
    end if
  end subroutine print_fit_summary

  subroutine write_simulation_csv(filename, simulation, status)
    character(len=*), intent(in) :: filename
    type(mfgarch_simulation), intent(in) :: simulation
    integer, intent(out) :: status
    integer :: unit, ios, i

    status = mfgarch_success
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      status = ios
      return
    end if
    write(unit,'(a)') 'day,return,covariate,low_frequency,tau,g,realized_variance,realized_variance_half_hour'
    do i = 1, size(simulation%returns)
      write(unit,'(i0,7(",",es24.16))') i, simulation%returns(i), simulation%covariate(i), &
        real(simulation%low_frequency_period(i),dp), simulation%tau(i), simulation%g(i), &
        simulation%realized_variance(i), simulation%realized_variance_half_hour(i)
    end do
    close(unit)
  end subroutine write_simulation_csv

end module mfgarch_io
