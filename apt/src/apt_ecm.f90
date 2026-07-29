! SPDX-License-Identifier: GPL-2.0-or-later
module apt_ecm
  use apt_kinds, only : dp
  use apt_regression, only : regression_result, hypothesis_result, &
    residual_diagnostics, fit_ols, linear_f_test, zero_coefficient_f_test, &
    compute_residual_diagnostics
  use apt_cointegration, only : apt_tar, apt_mtar
  implicit none
  private

  integer, parameter, public :: apt_linear = 0

  type, public :: ecm_fit_result
    integer :: lag = 1
    integer :: model = apt_linear
    integer :: n_total = 0
    integer :: n_regression = 0
    integer :: status = 0
    logical :: asymmetric = .false.
    logical :: split = .false.
    real(dp) :: threshold = 0.0_dp
    real(dp), allocatable :: long_run_coefficients(:)
    real(dp), allocatable :: long_run_residuals(:)
    real(dp), allocatable :: response_x(:)
    real(dp), allocatable :: response_y(:)
    real(dp), allocatable :: independent(:,:)
    type(regression_result) :: equation_x
    type(regression_result) :: equation_y
  end type ecm_fit_result

  type, public :: ecm_equation_diagnostics
    real(dp) :: r_squared = 0.0_dp
    real(dp) :: adjusted_r_squared = 0.0_dp
    real(dp) :: f_statistic = 0.0_dp
    real(dp) :: durbin_watson = 0.0_dp
    real(dp) :: durbin_watson_p = 1.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    real(dp) :: ljung_box_4_p = 1.0_dp
    real(dp) :: ljung_box_8_p = 1.0_dp
    real(dp) :: ljung_box_12_p = 1.0_dp
  end type ecm_equation_diagnostics

  type, public :: ecm_diagnostics_result
    integer :: status = 0
    type(ecm_equation_diagnostics) :: equation_x
    type(ecm_equation_diagnostics) :: equation_y
  end type ecm_diagnostics_result

  type, public :: paired_hypothesis_result
    character(len=96) :: description = ''
    type(hypothesis_result) :: equation_x
    type(hypothesis_result) :: equation_y
  end type paired_hypothesis_result

  type, public :: ecm_asymmetry_test_result
    integer :: status = 0
    integer :: ntests = 0
    type(paired_hypothesis_result), allocatable :: tests(:)
  end type ecm_asymmetry_test_result

  public :: ecm_symmetric_fit, ecm_asymmetric_fit, ecm_asymmetry_tests
  public :: ecm_diagnostics
  public :: ecmSymFit, ecmAsyFit, ecmAsyTest, ecmDiag

  interface ecmSymFit
    module procedure ecm_symmetric_fit
  end interface
  interface ecmAsyFit
    module procedure ecm_asymmetric_fit
  end interface
  interface ecmAsyTest
    module procedure ecm_asymmetry_tests
  end interface
  interface ecmDiag
    module procedure ecm_diagnostics
  end interface

contains

  subroutine fit_long_run(y, x, coefficients, residuals, status)
    real(dp), intent(in) :: y(:), x(:)
    real(dp), allocatable, intent(out) :: coefficients(:), residuals(:)
    integer, intent(out) :: status
    real(dp), allocatable :: design(:,:)
    type(regression_result) :: fit
    integer :: n
    n = size(y)
    if (size(x) /= n) then
      status = 1
      return
    end if
    allocate(design(n,2))
    design(:,1) = 1.0_dp
    design(:,2) = x
    call fit_ols(y, design, fit)
    if (fit%status /= 0) then
      status = 2
      return
    end if
    coefficients = fit%coefficients
    residuals = fit%residuals
    status = 0
  end subroutine fit_long_run

  subroutine ecm_symmetric_fit(y, x, result, lag)
    real(dp), intent(in) :: y(:), x(:)
    type(ecm_fit_result), intent(out) :: result
    integer, intent(in), optional :: lag
    integer :: n, l, t0, m, t, j, row, p, istat
    real(dp), allocatable :: design(:,:)

    n = size(y)
    l = 1
    if (present(lag)) l = lag
    result%lag = l
    result%n_total = n
    result%model = apt_linear
    result%asymmetric = .false.
    result%split = .false.
    if (size(x) /= n .or. l < 1 .or. n <= l + 2) then
      result%status = 1
      return
    end if
    call fit_long_run(y, x, result%long_run_coefficients, &
      result%long_run_residuals, istat)
    if (istat /= 0) then
      result%status = 2
      return
    end if
    t0 = l + 2
    m = n - t0 + 1
    p = 2 * l + 2
    result%n_regression = m
    allocate(result%response_x(m), result%response_y(m), &
      result%independent(m,p-1), design(m,p))
    do t = t0, n
      row = t - t0 + 1
      result%response_x(row) = x(t) - x(t-1)
      result%response_y(row) = y(t) - y(t-1)
      do j = 1, l
        result%independent(row,j) = x(t-j) - x(t-j-1)
        result%independent(row,l+j) = y(t-j) - y(t-j-1)
      end do
      result%independent(row,2*l+1) = result%long_run_residuals(t-1)
    end do
    design(:,1) = 1.0_dp
    design(:,2:p) = result%independent
    call fit_ols(result%response_x, design, result%equation_x)
    call fit_ols(result%response_y, design, result%equation_y)
    if (result%equation_x%status /= 0 .or. result%equation_y%status /= 0) then
      result%status = 3
      return
    end if
    result%status = 0
  end subroutine ecm_symmetric_fit

  subroutine ecm_asymmetric_fit(y, x, result, lag, split, model, threshold)
    real(dp), intent(in) :: y(:), x(:)
    type(ecm_fit_result), intent(out) :: result
    integer, intent(in), optional :: lag, model
    logical, intent(in), optional :: split
    real(dp), intent(in), optional :: threshold
    integer :: n, l, mdl, t0, m, t, j, row, q, p, istat, col
    logical :: spl
    real(dp) :: th, ect, indicator, dxlag, dylag
    real(dp), allocatable :: design(:,:)

    n = size(y)
    l = 1
    if (present(lag)) l = lag
    spl = .true.
    if (present(split)) spl = split
    mdl = apt_linear
    if (present(model)) mdl = model
    th = 0.0_dp
    if (present(threshold)) th = threshold
    result%lag = l
    result%split = spl
    result%model = mdl
    result%threshold = th
    result%asymmetric = .true.
    result%n_total = n
    if (size(x) /= n .or. l < 1 .or. n <= l + 2 .or. &
      (mdl /= apt_linear .and. mdl /= apt_tar .and. mdl /= apt_mtar)) then
      result%status = 1
      return
    end if
    call fit_long_run(y, x, result%long_run_coefficients, &
      result%long_run_residuals, istat)
    if (istat /= 0) then
      result%status = 2
      return
    end if
    t0 = l + 2
    if (mdl == apt_mtar) t0 = max(t0, 3)
    m = n - t0 + 1
    if (spl) then
      q = 4*l + 2
    else
      q = 2*l + 2
    end if
    p = q + 1
    result%n_regression = m
    allocate(result%response_x(m), result%response_y(m), &
      result%independent(m,q), design(m,p))
    result%independent = 0.0_dp

    do t = t0, n
      row = t - t0 + 1
      result%response_x(row) = x(t) - x(t-1)
      result%response_y(row) = y(t) - y(t-1)
      ect = result%long_run_residuals(t-1)
      select case (mdl)
      case (apt_linear)
        if (ect >= 0.0_dp) then
          indicator = 1.0_dp
        else
          indicator = 0.0_dp
        end if
      case (apt_tar)
        if (ect >= th) then
          indicator = 1.0_dp
        else
          indicator = 0.0_dp
        end if
      case (apt_mtar)
        if (result%long_run_residuals(t-1) - &
            result%long_run_residuals(t-2) >= th) then
          indicator = 1.0_dp
        else
          indicator = 0.0_dp
        end if
      end select
      if (spl) then
        do j = 1, l
          dxlag = x(t-j) - x(t-j-1)
          dylag = y(t-j) - y(t-j-1)
          result%independent(row,j) = max(dxlag, 0.0_dp)
          result%independent(row,l+j) = min(dxlag, 0.0_dp)
          result%independent(row,2*l+j) = max(dylag, 0.0_dp)
          result%independent(row,3*l+j) = min(dylag, 0.0_dp)
        end do
        result%independent(row,4*l+1) = ect * indicator
        result%independent(row,4*l+2) = ect * (1.0_dp-indicator)
      else
        do j = 1, l
          result%independent(row,j) = x(t-j) - x(t-j-1)
          result%independent(row,l+j) = y(t-j) - y(t-j-1)
        end do
        result%independent(row,2*l+1) = ect * indicator
        result%independent(row,2*l+2) = ect * (1.0_dp-indicator)
      end if
    end do
    design(:,1) = 1.0_dp
    do col = 1, q
      design(:,col+1) = result%independent(:,col)
    end do
    call fit_ols(result%response_x, design, result%equation_x)
    call fit_ols(result%response_y, design, result%equation_y)
    if (result%equation_x%status /= 0 .or. result%equation_y%status /= 0) then
      result%status = 3
      return
    end if
    result%status = 0
  end subroutine ecm_asymmetric_fit

  subroutine ecm_diagnostics(model, result)
    type(ecm_fit_result), intent(in) :: model
    type(ecm_diagnostics_result), intent(out) :: result
    type(residual_diagnostics) :: dx, dy
    if (model%status /= 0) then
      result%status = 1
      return
    end if
    call compute_residual_diagnostics(model%equation_x%residuals, dx)
    call compute_residual_diagnostics(model%equation_y%residuals, dy)
    call fill_equation_diag(model%equation_x, dx, result%equation_x)
    call fill_equation_diag(model%equation_y, dy, result%equation_y)
    result%status = 0
  end subroutine ecm_diagnostics

  subroutine fill_equation_diag(fit, rd, out)
    type(regression_result), intent(in) :: fit
    type(residual_diagnostics), intent(in) :: rd
    type(ecm_equation_diagnostics), intent(out) :: out
    out%r_squared = fit%r_squared
    out%adjusted_r_squared = fit%adjusted_r_squared
    out%f_statistic = fit%f_statistic
    out%durbin_watson = rd%durbin_watson
    out%durbin_watson_p = rd%durbin_watson_p
    out%aic = fit%aic
    out%bic = fit%bic
    out%ljung_box_4_p = rd%ljung_box_4_p
    out%ljung_box_8_p = rd%ljung_box_8_p
    out%ljung_box_12_p = rd%ljung_box_12_p
  end subroutine fill_equation_diag

  subroutine ecm_asymmetry_tests(model, result)
    type(ecm_fit_result), intent(in) :: model
    type(ecm_asymmetry_test_result), intent(out) :: result
    integer :: l, p, ntest, pos, j
    real(dp), allocatable :: r(:,:), rhs(:)
    integer, allocatable :: idx(:)

    if (model%status /= 0 .or. .not. model%asymmetric) then
      result%status = 1
      return
    end if
    l = model%lag
    p = model%equation_x%ncoef
    if (model%split) then
      ntest = 2*l + 5
    else
      ntest = 3
    end if
    result%ntests = ntest
    allocate(result%tests(ntest))

    allocate(r(1,p), rhs(1))
    r = 0.0_dp
    rhs = 0.0_dp
    r(1,p-1) = 1.0_dp
    r(1,p) = -1.0_dp
    result%tests(1)%description = 'Equilibrium adjustment path symmetry'
    call linear_f_test(model%equation_x, r, rhs, result%tests(1)%equation_x)
    call linear_f_test(model%equation_y, r, rhs, result%tests(1)%equation_y)
    deallocate(r, rhs)

    if (model%split) then
      allocate(idx(2*l))
      idx = [(j, j=2,2*l+1)]
    else
      allocate(idx(l))
      idx = [(j, j=2,l+1)]
    end if
    result%tests(2)%description = 'x does not Granger cause equation'
    call zero_coefficient_f_test(model%equation_x, idx, result%tests(2)%equation_x)
    call zero_coefficient_f_test(model%equation_y, idx, result%tests(2)%equation_y)
    if (model%split) then
      idx = [(j, j=2*l+2,4*l+1)]
    else
      idx = [(j, j=l+2,2*l+1)]
    end if
    result%tests(3)%description = 'y does not Granger cause equation'
    call zero_coefficient_f_test(model%equation_x, idx, result%tests(3)%equation_x)
    call zero_coefficient_f_test(model%equation_y, idx, result%tests(3)%equation_y)
    deallocate(idx)

    if (model%split) then
      pos = 3
      allocate(r(1,p), rhs(1))
      rhs = 0.0_dp
      do j = 1, l
        pos = pos + 1
        r = 0.0_dp
        r(1,1+j) = 1.0_dp
        r(1,1+l+j) = -1.0_dp
        write(result%tests(pos)%description,'(a,i0)') &
          'Distributed-lag symmetry for x at lag ', j
        call linear_f_test(model%equation_x, r, rhs, result%tests(pos)%equation_x)
        call linear_f_test(model%equation_y, r, rhs, result%tests(pos)%equation_y)
      end do
      do j = 1, l
        pos = pos + 1
        r = 0.0_dp
        r(1,1+2*l+j) = 1.0_dp
        r(1,1+3*l+j) = -1.0_dp
        write(result%tests(pos)%description,'(a,i0)') &
          'Distributed-lag symmetry for y at lag ', j
        call linear_f_test(model%equation_x, r, rhs, result%tests(pos)%equation_x)
        call linear_f_test(model%equation_y, r, rhs, result%tests(pos)%equation_y)
      end do
      pos = pos + 1
      r = 0.0_dp
      do j = 1, l
        r(1,1+j) = 1.0_dp
        r(1,1+l+j) = -1.0_dp
      end do
      result%tests(pos)%description = 'Cumulative x distributed-lag symmetry'
      call linear_f_test(model%equation_x, r, rhs, result%tests(pos)%equation_x)
      call linear_f_test(model%equation_y, r, rhs, result%tests(pos)%equation_y)
      pos = pos + 1
      r = 0.0_dp
      do j = 1, l
        r(1,1+2*l+j) = 1.0_dp
        r(1,1+3*l+j) = -1.0_dp
      end do
      result%tests(pos)%description = 'Cumulative y distributed-lag symmetry'
      call linear_f_test(model%equation_x, r, rhs, result%tests(pos)%equation_x)
      call linear_f_test(model%equation_y, r, rhs, result%tests(pos)%equation_y)
    end if
    result%status = 0
  end subroutine ecm_asymmetry_tests

end module apt_ecm
