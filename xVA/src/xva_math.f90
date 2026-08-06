module xva_math
  use trading, only : dp
  implicit none
  private

  public :: normal_pdf
  public :: inverse_normal_cdf
  public :: quadratic_form
  public :: same_text
  public :: uppercase

contains

  pure real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

    value = exp(-0.5_dp * x * x) / sqrt(two_pi)
  end function normal_pdf

  pure real(dp) function inverse_normal_cdf(probability) result(value)
    real(dp), intent(in) :: probability
    real(dp), parameter :: a1 = -3.969683028665376e+01_dp
    real(dp), parameter :: a2 =  2.209460984245205e+02_dp
    real(dp), parameter :: a3 = -2.759285104469687e+02_dp
    real(dp), parameter :: a4 =  1.383577518672690e+02_dp
    real(dp), parameter :: a5 = -3.066479806614716e+01_dp
    real(dp), parameter :: a6 =  2.506628277459239e+00_dp
    real(dp), parameter :: b1 = -5.447609879822406e+01_dp
    real(dp), parameter :: b2 =  1.615858368580409e+02_dp
    real(dp), parameter :: b3 = -1.556989798598866e+02_dp
    real(dp), parameter :: b4 =  6.680131188771972e+01_dp
    real(dp), parameter :: b5 = -1.328068155288572e+01_dp
    real(dp), parameter :: c1 = -7.784894002430293e-03_dp
    real(dp), parameter :: c2 = -3.223964580411365e-01_dp
    real(dp), parameter :: c3 = -2.400758277161838e+00_dp
    real(dp), parameter :: c4 = -2.549732539343734e+00_dp
    real(dp), parameter :: c5 =  4.374664141464968e+00_dp
    real(dp), parameter :: c6 =  2.938163982698783e+00_dp
    real(dp), parameter :: d1 =  7.784695709041462e-03_dp
    real(dp), parameter :: d2 =  3.224671290700398e-01_dp
    real(dp), parameter :: d3 =  2.445134137142996e+00_dp
    real(dp), parameter :: d4 =  3.754408661907416e+00_dp
    real(dp), parameter :: p_low = 0.02425_dp
    real(dp), parameter :: p_high = 1.0_dp - p_low
    real(dp) :: error_value
    real(dp) :: q
    real(dp) :: r

    if (probability <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    else if (probability >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    end if

    if (probability < p_low) then
      q = sqrt(-2.0_dp * log(probability))
      value = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
    else if (probability <= p_high) then
      q = probability - 0.5_dp
      r = q * q
      value = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
        (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - probability))
      value = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
        ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
    end if

    error_value = 0.5_dp * erfc(-value / sqrt(2.0_dp)) - probability
    value = value - error_value / normal_pdf(value)
  end function inverse_normal_cdf

  pure real(dp) function quadratic_form(vector, matrix) result(value)
    real(dp), intent(in) :: vector(:)
    real(dp), intent(in) :: matrix(:,:)
    real(dp) :: temp(size(vector))

    if (size(matrix, 1) /= size(vector) .or. &
        size(matrix, 2) /= size(vector)) then
      value = huge(1.0_dp)
      return
    end if
    temp = matmul(matrix, vector)
    value = dot_product(vector, temp)
  end function quadratic_form

  pure logical function same_text(left, right) result(equal)
    character(len=*), intent(in) :: left
    character(len=*), intent(in) :: right

    equal = uppercase(trim(adjustl(left))) == uppercase(trim(adjustl(right)))
  end function same_text

  pure function uppercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: code
    integer :: i

    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        value(i:i) = achar(code - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

end module xva_math
