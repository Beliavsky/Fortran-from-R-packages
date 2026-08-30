! SPDX-License-Identifier: GPL-2.0-only
module multcomp_math
  use multcomp_kinds, only : dp
  use mvtnorm, only : jacobi_eigen, normal_cdf, student_t_cdf, &
    regularized_beta, chi_square_cdf
  implicit none
  private

  public :: covariance_to_correlation_mc
  public :: generalized_inverse
  public :: p_adjust
  public :: normal_pvalue
  public :: student_pvalue
  public :: chi_square_upper
  public :: f_upper
  public :: stable_order

contains

  subroutine covariance_to_correlation_mc(covariance, correlation, sd, ok)
    real(dp), intent(in) :: covariance(:, :) !! Symmetric covariance matrix to standardize.
    real(dp), allocatable, intent(out) :: correlation(:, :) !! Correlation matrix with unit diagonal.
    real(dp), allocatable, intent(out) :: sd(:) !! Marginal standard deviations from the covariance diagonal.
    logical, intent(out) :: ok !! True when every covariance diagonal entry is positive.

    integer :: i
    integer :: j
    integer :: n

    n = size(covariance, 1)
    allocate(correlation(n, n), sd(n))
    correlation = 0.0_dp
    sd = 0.0_dp
    ok = size(covariance, 2) == n
    if (.not. ok) return

    do i = 1, n
      if (covariance(i, i) <= 0.0_dp) then
        ok = .false.
        return
      end if
      sd(i) = sqrt(covariance(i, i))
    end do

    do i = 1, n
      do j = 1, n
        correlation(i, j) = covariance(i, j) / (sd(i) * sd(j))
      end do
      correlation(i, i) = 1.0_dp
    end do
  end subroutine covariance_to_correlation_mc

  subroutine generalized_inverse(a, ainv, rank, ok, tolerance)
    real(dp), intent(in) :: a(:, :) !! Matrix whose Moore-Penrose generalized inverse is required.
    real(dp), allocatable, intent(out) :: ainv(:, :) !! Generalized inverse with shape (ncol(a), nrow(a)).
    integer, intent(out) :: rank !! Numerical rank under the MASS::ginv-style tolerance rule.
    logical, intent(out) :: ok !! True when the symmetric eigensolver converged.
    real(dp), intent(in), optional :: tolerance !! Relative singular-value cutoff; default sqrt(machine epsilon).

    real(dp), allocatable :: ata(:, :)
    real(dp), allocatable :: av(:)
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: vectors(:, :)
    real(dp) :: cutoff
    real(dp) :: lambda_max
    real(dp) :: tol
    integer :: i
    integer :: m
    integer :: n

    m = size(a, 1)
    n = size(a, 2)
    allocate(ainv(n, m))
    ainv = 0.0_dp
    rank = 0
    ok = size(a, 1) >= 1 .and. size(a, 2) >= 1
    if (.not. ok) return

    ata = matmul(transpose(a), a)
    call jacobi_eigen(ata, values, vectors, ok)
    if (.not. ok) return

    lambda_max = max(0.0_dp, maxval(values))
    tol = sqrt(epsilon(1.0_dp))
    if (present(tolerance)) tol = max(0.0_dp, tolerance)
    cutoff = tol * tol * lambda_max
    allocate(av(m))

    do i = 1, n
      if (values(i) > cutoff) then
        av = matmul(a, vectors(:, i))
        ainv = ainv + spread(vectors(:, i), 2, m) * &
          spread(av / values(i), 1, n)
        rank = rank + 1
      end if
    end do
  end subroutine generalized_inverse

  subroutine stable_order(x, order, decreasing)
    real(dp), intent(in) :: x(:) !! Values to order, with original index breaking exact ties.
    integer, allocatable, intent(out) :: order(:) !! One-based permutation of x in requested order.
    logical, intent(in), optional :: decreasing !! If true, order largest to smallest; default false.

    logical :: desc
    integer :: i
    integer :: j
    integer :: key
    integer :: n

    n = size(x)
    allocate(order(n))
    order = [(i, i = 1, n)]
    desc = .false.
    if (present(decreasing)) desc = decreasing

    do i = 2, n
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (.not. comes_after(order(j), key)) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do

  contains

    logical function comes_after(left_index, right_index) result(after)
      integer, intent(in) :: left_index !! Current permutation index being compared.
      integer, intent(in) :: right_index !! Candidate permutation index being inserted.

      if (desc) then
        if (x(left_index) < x(right_index)) then
          after = .true.
        else if (x(left_index) > x(right_index)) then
          after = .false.
        else
          after = left_index > right_index
        end if
      else
        if (x(left_index) > x(right_index)) then
          after = .true.
        else if (x(left_index) < x(right_index)) then
          after = .false.
        else
          after = left_index > right_index
        end if
      end if
    end function comes_after

  end subroutine stable_order

  recursive subroutine p_adjust(p, method, adjusted, ok)
    real(dp), intent(in) :: p(:) !! Raw p-values, normally constrained to the interval [0,1].
    character(len=*), intent(in) :: method !! Adjustment name: none, bonferroni, holm, hochberg, hommel, BH/fdr, or BY.
    real(dp), allocatable, intent(out) :: adjusted(:) !! Multiplicity-adjusted p-values in original order.
    logical, intent(out) :: ok !! True when the requested method is recognized.

    character(len=:), allocatable :: key
    integer, allocatable :: ord(:)
    integer :: i
    integer :: j
    integer :: m
    real(dp), allocatable :: ps(:)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: pa(:)
    real(dp) :: harmonic
    real(dp) :: q1

    m = size(p)
    allocate(adjusted(m))
    if (m == 0) then
      ok = .true.
      return
    end if

    key = lowercase(trim(adjustl(method)))
    ok = .true.

    select case (key)
    case ('none')
      adjusted = min(1.0_dp, max(0.0_dp, p))
    case ('bonferroni')
      adjusted = min(1.0_dp, real(m, dp) * p)
    case ('holm')
      call stable_order(p, ord)
      allocate(ps(m), q(m))
      ps = p(ord)
      q(1) = real(m, dp) * ps(1)
      do i = 2, m
        q(i) = max(q(i - 1), real(m - i + 1, dp) * ps(i))
      end do
      do i = 1, m
        adjusted(ord(i)) = min(1.0_dp, q(i))
      end do
    case ('hochberg')
      call stable_order(p, ord, decreasing=.true.)
      allocate(ps(m), q(m))
      ps = p(ord)
      q(1) = ps(1)
      do i = 2, m
        q(i) = min(q(i - 1), real(i, dp) * ps(i))
      end do
      do i = 1, m
        adjusted(ord(i)) = min(1.0_dp, q(i))
      end do
    case ('bh', 'fdr')
      call stable_order(p, ord, decreasing=.true.)
      allocate(ps(m), q(m))
      ps = p(ord)
      q(1) = ps(1)
      do i = 2, m
        q(i) = min(q(i - 1), real(m, dp) / real(m - i + 1, dp) * ps(i))
      end do
      do i = 1, m
        adjusted(ord(i)) = min(1.0_dp, q(i))
      end do
    case ('by')
      harmonic = 0.0_dp
      do i = 1, m
        harmonic = harmonic + 1.0_dp / real(i, dp)
      end do
      call stable_order(p, ord, decreasing=.true.)
      allocate(ps(m), q(m))
      ps = p(ord)
      q(1) = harmonic * ps(1)
      do i = 2, m
        q(i) = min(q(i - 1), harmonic * real(m, dp) / &
          real(m - i + 1, dp) * ps(i))
      end do
      do i = 1, m
        adjusted(ord(i)) = min(1.0_dp, q(i))
      end do
    case ('hommel')
      if (m <= 2) then
        call p_adjust(p, 'holm', adjusted, ok)
        return
      end if
      call stable_order(p, ord)
      allocate(ps(m), q(m), pa(m))
      ps = p(ord)
      q1 = huge(1.0_dp)
      do i = 1, m
        q1 = min(q1, real(m, dp) * ps(i) / real(i, dp))
      end do
      q = q1
      pa = q1
      do j = m - 1, 2, -1
        q1 = huge(1.0_dp)
        do i = m - j + 2, m
          q1 = min(q1, real(j, dp) * ps(i) / real(i - (m - j), dp))
        end do
        do i = 1, m - j + 1
          q(i) = min(real(j, dp) * ps(i), q1)
        end do
        do i = m - j + 2, m
          q(i) = q(m - j + 1)
        end do
        pa = max(pa, q)
      end do
      pa = max(pa, ps)
      do i = 1, m
        adjusted(ord(i)) = min(1.0_dp, pa(i))
      end do
    case default
      adjusted = p
      ok = .false.
    end select
  end subroutine p_adjust

  real(dp) function normal_pvalue(statistic, alternative) result(p)
    real(dp), intent(in) :: statistic !! Standard-normal test statistic.
    integer, intent(in) :: alternative !! Alternative code: 1 two-sided, 2 less, 3 greater.

    select case (alternative)
    case (2)
      p = normal_cdf(statistic)
    case (3)
      p = 1.0_dp - normal_cdf(statistic)
    case default
      p = 2.0_dp * (1.0_dp - normal_cdf(abs(statistic)))
    end select
    p = min(1.0_dp, max(0.0_dp, p))
  end function normal_pvalue

  real(dp) function student_pvalue(statistic, df, alternative) result(p)
    real(dp), intent(in) :: statistic !! Student t test statistic.
    real(dp), intent(in) :: df !! Positive Student t degrees of freedom.
    integer, intent(in) :: alternative !! Alternative code: 1 two-sided, 2 less, 3 greater.

    select case (alternative)
    case (2)
      p = student_t_cdf(statistic, df)
    case (3)
      p = 1.0_dp - student_t_cdf(statistic, df)
    case default
      p = 2.0_dp * (1.0_dp - student_t_cdf(abs(statistic), df))
    end select
    p = min(1.0_dp, max(0.0_dp, p))
  end function student_pvalue

  real(dp) function chi_square_upper(statistic, df) result(p)
    real(dp), intent(in) :: statistic !! Nonnegative chi-square statistic.
    real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.

    p = 1.0_dp - chi_square_cdf(max(0.0_dp, statistic), df)
    p = min(1.0_dp, max(0.0_dp, p))
  end function chi_square_upper

  real(dp) function f_upper(statistic, df1, df2) result(p)
    real(dp), intent(in) :: statistic !! Nonnegative F statistic.
    real(dp), intent(in) :: df1 !! Positive numerator degrees of freedom.
    real(dp), intent(in) :: df2 !! Positive denominator degrees of freedom.

    real(dp) :: z

    if (statistic <= 0.0_dp) then
      p = 1.0_dp
      return
    end if
    z = df1 * statistic / (df1 * statistic + df2)
    p = 1.0_dp - regularized_beta(z, 0.5_dp * df1, 0.5_dp * df2)
    p = min(1.0_dp, max(0.0_dp, p))
  end function f_upper

  pure function lowercase(text) result(out)
    character(len=*), intent(in) :: text !! Text to normalize for case-insensitive method matching.
    character(len=len(text)) :: out

    integer :: c
    integer :: i

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function lowercase

end module multcomp_math
