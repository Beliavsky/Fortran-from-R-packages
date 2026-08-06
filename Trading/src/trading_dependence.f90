module trading_dependence
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use trading_kinds, only : dp
  use trading_stats, only : correlation_dp, correlation_matrix, sd_dp, variance_dp
  implicit none
  private

  public :: angular_distance
  public :: chebyshev_distance
  public :: sample_entropy
  public :: cross_sample_entropy
  public :: normalized_cross_sample_entropy
  public :: variation_of_information
  public :: information_adjusted_correlation
  public :: information_adjusted_beta

contains

  subroutine angular_distance(returns_matrix, distance, long_short)
    real(dp), intent(in) :: returns_matrix(:, :)
    real(dp), intent(out) :: distance(:, :)
    logical, intent(in), optional :: long_short
    real(dp), allocatable :: corr(:, :)
    logical :: ignore_sign
    integer :: p

    p = size(returns_matrix, 2)
    if (p < 2) error stop "angular_distance: at least two columns are required"
    if (size(distance, 1) /= p .or. size(distance, 2) /= p) then
      error stop "angular_distance: result has the wrong shape"
    end if

    ignore_sign = .false.
    if (present(long_short)) ignore_sign = long_short

    allocate(corr(p, p))
    call correlation_matrix(returns_matrix, corr)
    if (ignore_sign) corr = abs(corr)
    distance = sqrt(max((1.0_dp - corr) / 2.0_dp, 0.0_dp))
  end subroutine angular_distance

  pure real(dp) function chebyshev_distance(x, y) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer :: n

    n = min(size(x), size(y))
    if (n == 0) then
      value = 0.0_dp
    else
      value = maxval(abs(x(:n) - y(:n)))
    end if
  end function chebyshev_distance

  real(dp) function sample_entropy(returns, m, tolerance) result(value)
    real(dp), intent(in) :: returns(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: z(:)
    real(dp) :: r
    real(dp) :: scale
    integer :: a_count
    integer :: b_count
    integer :: embed
    integer :: i
    integer :: j
    integer :: n

    n = size(returns)
    embed = 2
    if (present(m)) embed = m
    r = 0.2_dp
    if (present(tolerance)) r = tolerance

    if (embed < 1) error stop "sample_entropy: m must be positive"
    if (n <= embed + 1) error stop "sample_entropy: input is too short"
    if (r <= 0.0_dp) error stop "sample_entropy: tolerance must be positive"

    scale = sd_dp(returns)
    if (scale <= tiny(1.0_dp)) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    allocate(z(n))
    z = (returns - sum(returns) / real(n, dp)) / scale

    b_count = 0
    do i = 1, n - embed
      do j = i + 1, n - embed + 1
        if (maxval(abs(z(i:i + embed - 1) - z(j:j + embed - 1))) < r) then
          b_count = b_count + 1
        end if
      end do
    end do

    a_count = 0
    do i = 1, n - embed - 1
      do j = i + 1, n - embed
        if (maxval(abs(z(i:i + embed) - z(j:j + embed))) < r) then
          a_count = a_count + 1
        end if
      end do
    end do

    if (a_count == 0 .or. b_count == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = -log(real(a_count, dp) / real(b_count, dp))
    end if
  end function sample_entropy

  real(dp) function cross_sample_entropy(x, y, m, tolerance) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: zx(:)
    real(dp), allocatable :: zy(:)
    real(dp) :: a_probability
    real(dp) :: b_probability
    real(dp) :: r
    real(dp) :: sx
    real(dp) :: sy
    integer :: a_count
    integer :: b_count
    integer :: embed
    integer :: i
    integer :: j
    integer :: n
    integer :: n_a
    integer :: n_b

    n = min(size(x), size(y))
    embed = 2
    if (present(m)) embed = m
    r = 0.2_dp
    if (present(tolerance)) r = tolerance

    if (embed < 1) error stop "cross_sample_entropy: m must be positive"
    if (n <= embed + 1) error stop "cross_sample_entropy: input is too short"
    if (r <= 0.0_dp) error stop "cross_sample_entropy: tolerance must be positive"

    sx = sd_dp(x(:n))
    sy = sd_dp(y(:n))
    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    allocate(zx(n), zy(n))
    zx = (x(:n) - sum(x(:n)) / real(n, dp)) / sx
    zy = (y(:n) - sum(y(:n)) / real(n, dp)) / sy

    n_b = n - embed + 1
    b_count = 0
    do i = 1, n_b
      do j = 1, n_b
        if (maxval(abs(zx(i:i + embed - 1) - zy(j:j + embed - 1))) < r) then
          b_count = b_count + 1
        end if
      end do
    end do
    b_probability = real(b_count, dp) / real(n_b * n_b, dp)

    n_a = n - embed
    a_count = 0
    do i = 1, n_a
      do j = 1, n_a
        if (maxval(abs(zx(i:i + embed) - zy(j:j + embed))) < r) then
          a_count = a_count + 1
        end if
      end do
    end do
    a_probability = real(a_count, dp) / real(n_a * n_a, dp)

    if (a_probability <= 0.0_dp .or. b_probability <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = -log(a_probability / b_probability)
    end if
  end function cross_sample_entropy

  real(dp) function normalized_cross_sample_entropy(x, y, m, tolerance) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    real(dp) :: cross_entropy
    real(dp) :: entropy_x
    real(dp) :: entropy_y

    cross_entropy = cross_sample_entropy(x, y, m, tolerance)
    entropy_x = sample_entropy(x, m, tolerance)
    entropy_y = sample_entropy(y, m, tolerance)
    value = cross_entropy / min(entropy_x, entropy_y)
  end function normalized_cross_sample_entropy

  real(dp) function variation_of_information(x, y, m, tolerance, normalized) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: normalized
    logical :: do_normalize
    real(dp) :: cross_entropy
    real(dp) :: entropy_x
    real(dp) :: entropy_y
    real(dp) :: joint_distance

    do_normalize = .true.
    if (present(normalized)) do_normalize = normalized

    cross_entropy = cross_sample_entropy(x, y, m, tolerance)
    entropy_x = sample_entropy(x, m, tolerance)
    entropy_y = sample_entropy(y, m, tolerance)

    value = entropy_x + entropy_y - 2.0_dp * cross_entropy
    if (do_normalize) then
      joint_distance = entropy_x + entropy_y - cross_entropy
      if (abs(joint_distance) <= tiny(1.0_dp)) then
        value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
        value = value / joint_distance
      end if
    end if
  end function variation_of_information

  real(dp) function information_adjusted_correlation(x, y, m, tolerance) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    real(dp) :: mutual_information_measure
    real(dp) :: pearson
    real(dp) :: magnitude

    mutual_information_measure = normalized_cross_sample_entropy(x, y, m, tolerance)
    pearson = correlation_dp(x, y)

    if (mutual_information_measure <= 0.0_dp) then
      value = 0.0_dp
      return
    end if

    magnitude = sqrt(max(1.0_dp - 2.0_dp**(-2.0_dp * mutual_information_measure), 0.0_dp))
    value = sign(magnitude, pearson)
  end function information_adjusted_correlation

  real(dp) function information_adjusted_beta(x, y, m, tolerance) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    integer, intent(in), optional :: m
    real(dp), intent(in), optional :: tolerance
    real(dp) :: var_y

    var_y = variance_dp(y)
    if (var_y <= tiny(1.0_dp)) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = information_adjusted_correlation(x, y, m, tolerance) * &
        sqrt(variance_dp(x) / var_y)
    end if
  end function information_adjusted_beta

end module trading_dependence
