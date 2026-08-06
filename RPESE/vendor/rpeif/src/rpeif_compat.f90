! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_compat
  use rpeif_kinds, only : dp
  use rpeif_types, only : rpeif_options
  use rpeif_influence, only : influence_from_data
  implicit none
  private
  public :: if_mean, if_sd, if_semisd, if_var, if_es, if_sr, if_sortino
  public :: if_downside_sharpe, if_es_ratio, if_var_ratio, if_rachev_ratio
  public :: if_robust_mean, if_lpm, if_omega_ratio
contains
  subroutine if_mean(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('mean', x, returns, values, options, status)
  end subroutine if_mean

  subroutine if_sd(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('sd', x, returns, values, options, status)
  end subroutine if_sd

  subroutine if_semisd(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('semisd', x, returns, values, options, status)
  end subroutine if_semisd

  subroutine if_var(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('var', x, returns, values, options, status)
  end subroutine if_var

  subroutine if_es(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('es', x, returns, values, options, status)
  end subroutine if_es

  subroutine if_sr(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('sr', x, returns, values, options, status)
  end subroutine if_sr

  subroutine if_sortino(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('sor', x, returns, values, options, status)
  end subroutine if_sortino

  subroutine if_downside_sharpe(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('dsr', x, returns, values, options, status)
  end subroutine if_downside_sharpe

  subroutine if_es_ratio(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('esratio', x, returns, values, options, status)
  end subroutine if_es_ratio

  subroutine if_var_ratio(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('varratio', x, returns, values, options, status)
  end subroutine if_var_ratio

  subroutine if_rachev_ratio(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('rachevratio', x, returns, values, options, status)
  end subroutine if_rachev_ratio

  subroutine if_robust_mean(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('robmean', x, returns, values, options, status)
  end subroutine if_robust_mean

  subroutine if_lpm(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('lpm', x, returns, values, options, status)
  end subroutine if_lpm

  subroutine if_omega_ratio(x, returns, values, options, status)
    real(dp), intent(in) :: x(:), returns(:)
    real(dp), allocatable, intent(out) :: values(:)
    type(rpeif_options), intent(in), optional :: options
    integer, intent(out), optional :: status
    call influence_from_data('omegaratio', x, returns, values, options, status)
  end subroutine if_omega_ratio
end module rpeif_compat
