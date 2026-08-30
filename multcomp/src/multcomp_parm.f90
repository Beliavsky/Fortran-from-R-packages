! SPDX-License-Identifier: GPL-2.0-only
module multcomp_parm
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type
  implicit none
  private

  public :: make_parm
  public :: parm_from_iid
  public :: block_diag_covariance

contains

  subroutine make_parm(coef, vcov, result, df)
    real(dp), intent(in) :: coef(:) !! Parameter estimates in the same order as vcov rows and columns.
    real(dp), intent(in) :: vcov(:, :) !! Symmetric covariance matrix of the parameter estimates.
    type(parm_type), intent(out) :: result !! Validated parameter object used by the GLHT routines.
    real(dp), intent(in), optional :: df !! Residual degrees of freedom; zero requests asymptotic normal inference.

    real(dp) :: scale
    real(dp) :: tolerance
    integer :: n

    n = size(coef)
    result%ok = .false.
    result%message = ''
    if (size(vcov, 1) /= n .or. size(vcov, 2) /= n) then
      result%message = 'coefficient and covariance dimensions do not match'
      return
    end if

    scale = max(1.0_dp, maxval(abs(vcov)))
    tolerance = sqrt(epsilon(1.0_dp)) * scale
    if (maxval(abs(vcov - transpose(vcov))) > tolerance) then
      result%message = 'covariance matrix is not symmetric within tolerance'
      return
    end if

    allocate(result%coef(n), result%vcov(n, n))
    result%coef = coef
    result%vcov = 0.5_dp * (vcov + transpose(vcov))
    result%df = 0.0_dp
    if (present(df)) result%df = df
    if (result%df < 0.0_dp) then
      result%message = 'degrees of freedom must be nonnegative'
      return
    end if

    result%ok = .true.
  end subroutine make_parm

  subroutine parm_from_iid(coef, iid, result, model_standard_error, df)
    real(dp), intent(in) :: coef(:) !! Parameter estimates corresponding to the rows of iid.
    real(dp), intent(in) :: iid(:, :) !! IID influence contributions, with parameters by observations.
    type(parm_type), intent(out) :: result !! Parameter object with covariance estimated from IID contributions.
    real(dp), intent(in), optional :: model_standard_error(:) !! Marginal model-based standard errors replacing IID scales.
    real(dp), intent(in), optional :: df !! Residual degrees of freedom; zero requests asymptotic normal inference.

    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: scale(:)
    integer :: i
    integer :: j
    integer :: nobs
    integer :: p

    p = size(coef)
    nobs = size(iid, 2)
    result%ok = .false.
    result%message = ''
    if (size(iid, 1) /= p .or. nobs < 1) then
      result%message = 'iid must have one row per coefficient and at least one observation'
      return
    end if

    covariance = matmul(iid, transpose(iid)) / real(nobs * nobs, dp)
    if (present(model_standard_error)) then
      if (size(model_standard_error) /= p) then
        result%message = 'model_standard_error length does not match coefficient count'
        return
      end if
      allocate(scale(p), correlation(p, p))
      do i = 1, p
        if (covariance(i, i) <= 0.0_dp) then
          result%message = 'IID covariance has a nonpositive diagonal entry'
          return
        end if
        scale(i) = sqrt(covariance(i, i))
      end do
      do i = 1, p
        do j = 1, p
          correlation(i, j) = covariance(i, j) / (scale(i) * scale(j))
        end do
      end do
      do i = 1, p
        do j = 1, p
          covariance(i, j) = correlation(i, j) * model_standard_error(i) * &
            model_standard_error(j)
        end do
      end do
    end if

    call make_parm(coef, covariance, result, df)
  end subroutine parm_from_iid

  subroutine block_diag_covariance(blocks, result, ok)
    type(parm_type), intent(in) :: blocks(:) !! Parameter blocks to concatenate assuming zero cross-block covariance.
    type(parm_type), intent(out) :: result !! Combined parameter object with block-diagonal covariance.
    logical, intent(out) :: ok !! True when all blocks are valid and concatenation succeeds.

    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: covariance(:, :)
    integer :: first
    integer :: i
    integer :: last
    integer :: total

    total = 0
    ok = .true.
    do i = 1, size(blocks)
      if (.not. blocks(i)%ok) ok = .false.
      total = total + size(blocks(i)%coef)
    end do
    if (.not. ok .or. total == 0) then
      result%ok = .false.
      result%message = 'all parameter blocks must be valid and nonempty'
      return
    end if

    allocate(coef(total), covariance(total, total))
    covariance = 0.0_dp
    first = 1
    do i = 1, size(blocks)
      last = first + size(blocks(i)%coef) - 1
      coef(first:last) = blocks(i)%coef
      covariance(first:last, first:last) = blocks(i)%vcov
      first = last + 1
    end do
    call make_parm(coef, covariance, result)
  end subroutine block_diag_covariance

end module multcomp_parm
