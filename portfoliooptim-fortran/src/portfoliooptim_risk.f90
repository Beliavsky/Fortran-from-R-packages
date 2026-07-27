! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_risk
  use portfoliooptim_kinds, only : dp
  use portfoliooptim_types, only : risk_result, risk_cvar, risk_dcvar, risk_lsad, risk_mad
  implicit none
  private
  public :: risk_post, risk_measure, risk_code

contains

  function risk_post(losses, probabilities, alpha) result(res)
    real(dp), intent(in) :: losses(:)
    real(dp), intent(in) :: probabilities(:)
    real(dp), intent(in) :: alpha
    type(risk_result) :: res
    real(dp), allocatable :: x(:), p(:)
    real(dp) :: psum, cumulative, tail
    integer :: n, i, index

    n = size(losses)
    if (n == 0 .or. size(probabilities) /= n) then
      res%message = 'loss and probability dimensions are inconsistent'
      return
    end if
    if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      res%message = 'alpha must lie strictly between zero and one'
      return
    end if
    if (any(probabilities < 0.0_dp)) then
      res%message = 'probabilities must be nonnegative'
      return
    end if
    psum = sum(probabilities)
    if (psum <= 0.0_dp) then
      res%message = 'probabilities must have positive sum'
      return
    end if

    allocate(x(n), p(n))
    x = losses
    p = probabilities / psum
    call sort_pairs(x, p)
    res%mean = dot_product(x, p)
    res%mad = dot_product(abs(x - res%mean), p)

    cumulative = 0.0_dp
    index = n
    do i = 1, n
      cumulative = cumulative + p(i)
      if (cumulative >= alpha) then
        index = i
        exit
      end if
    end do
    res%var = x(index)
    tail = x(index) * max(0.0_dp, cumulative - alpha)
    if (index < n) tail = tail + dot_product(x(index + 1:n), p(index + 1:n))
    res%cvar = tail / (1.0_dp - alpha)
    res%ok = .true.
    res%message = 'ok'
  end function risk_post

  function risk_measure(losses, probabilities, alpha, code, portfolio_mean) result(value)
    real(dp), intent(in) :: losses(:), probabilities(:), alpha
    integer, intent(in) :: code
    real(dp), intent(in), optional :: portfolio_mean
    real(dp) :: value
    type(risk_result) :: summary
    real(dp) :: mean_value

    summary = risk_post(losses, probabilities, alpha)
    if (.not. summary%ok) then
      value = huge(1.0_dp)
      return
    end if
    mean_value = -summary%mean
    if (present(portfolio_mean)) mean_value = portfolio_mean
    select case (code)
    case (risk_cvar)
      value = summary%cvar
    case (risk_dcvar)
      value = summary%cvar + mean_value
    case (risk_lsad)
      value = 0.5_dp * summary%mad
    case (risk_mad)
      value = summary%mad
    case default
      value = huge(1.0_dp)
    end select
  end function risk_measure

  integer function risk_code(name) result(code)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: upper

    upper = uppercase(trim(adjustl(name)))
    select case (upper)
    case ('CVAR')
      code = risk_cvar
    case ('DCVAR')
      code = risk_dcvar
    case ('LSAD')
      code = risk_lsad
    case ('MAD')
      code = risk_mad
    case default
      code = 0
    end select
  end function risk_code

  subroutine sort_pairs(x, p)
    real(dp), intent(inout) :: x(:), p(:)
    integer :: i, j
    real(dp) :: key_x, key_p

    do i = 2, size(x)
      key_x = x(i)
      key_p = p(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key_x) exit
        x(j + 1) = x(j)
        p(j + 1) = p(j)
        j = j - 1
      end do
      x(j + 1) = key_x
      p(j + 1) = key_p
    end do
  end subroutine sort_pairs

  pure function uppercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
    end do
  end function uppercase

end module portfoliooptim_risk
