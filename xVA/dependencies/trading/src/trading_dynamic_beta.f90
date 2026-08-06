module trading_dynamic_beta
  use trading_kinds, only : dp
  implicit none
  private

  type, public :: dynamic_beta_result
    real(dp), allocatable :: filtered_state(:, :)
    real(dp), allocatable :: filtered_covariance(:, :, :)
    real(dp), allocatable :: smoothed_state(:, :)
    real(dp), allocatable :: smoothed_covariance(:, :, :)
    real(dp) :: transition_covariance(2, 2) = 0.0_dp
    real(dp) :: observation_variance = 0.0_dp
    real(dp) :: initial_covariance(2, 2) = 0.0_dp
  end type dynamic_beta_result

  public :: dynamic_beta

contains

  subroutine dynamic_beta(fund_returns, benchmark_returns, result, em_iterations, scale_percent)
    real(dp), intent(in) :: fund_returns(:)
    real(dp), intent(in) :: benchmark_returns(:)
    type(dynamic_beta_result), intent(out) :: result
    integer, intent(in), optional :: em_iterations
    logical, intent(in), optional :: scale_percent
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: predicted_mean(:, :)
    real(dp), allocatable :: predicted_covariance(:, :, :)
    real(dp), allocatable :: filtered_mean(:, :)
    real(dp), allocatable :: filtered_covariance(:, :, :)
    real(dp), allocatable :: smoothed_mean(:, :)
    real(dp), allocatable :: smoothed_covariance(:, :, :)
    real(dp), allocatable :: smoother_gain(:, :, :)
    real(dp) :: initial_mean(2)
    real(dp) :: initial_covariance(2, 2)
    real(dp) :: transition_covariance(2, 2)
    real(dp) :: observation_variance
    real(dp) :: delta
    logical :: do_scale
    integer :: iteration
    integer :: iterations
    integer :: n

    n = min(size(fund_returns), size(benchmark_returns))
    if (n < 3) error stop "dynamic_beta: at least three observations are required"

    iterations = 100
    if (present(em_iterations)) iterations = max(0, em_iterations)
    do_scale = .true.
    if (present(scale_percent)) do_scale = scale_percent

    allocate(y(n), x(n))
    y = fund_returns(:n)
    x = benchmark_returns(:n)
    if (do_scale) then
      y = 100.0_dp * y
      x = 100.0_dp * x
    end if

    allocate(predicted_mean(2, n), predicted_covariance(2, 2, n))
    allocate(filtered_mean(2, n), filtered_covariance(2, 2, n))
    allocate(smoothed_mean(2, n), smoothed_covariance(2, 2, n))
    allocate(smoother_gain(2, 2, max(n - 1, 1)))

    delta = 1.0e-5_dp
    initial_mean = 0.0_dp
    initial_covariance = 1.0_dp
    initial_covariance(1, 1) = initial_covariance(1, 1) + 1.0e-8_dp
    initial_covariance(2, 2) = initial_covariance(2, 2) + 1.0e-8_dp
    transition_covariance = 0.0_dp
    transition_covariance(1, 1) = delta / (1.0_dp - delta)
    transition_covariance(2, 2) = delta / (1.0_dp - delta)
    observation_variance = 1.0_dp

    do iteration = 1, iterations
      call filter_and_smooth(y, x, initial_mean, initial_covariance, &
        transition_covariance, observation_variance, predicted_mean, &
        predicted_covariance, filtered_mean, filtered_covariance, &
        smoothed_mean, smoothed_covariance, smoother_gain)
      call em_update(y, x, initial_mean, initial_covariance, &
        transition_covariance, observation_variance, smoothed_mean, &
        smoothed_covariance, smoother_gain)
    end do

    call filter_and_smooth(y, x, initial_mean, initial_covariance, &
      transition_covariance, observation_variance, predicted_mean, &
      predicted_covariance, filtered_mean, filtered_covariance, &
      smoothed_mean, smoothed_covariance, smoother_gain)

    allocate(result%filtered_state(2, n))
    allocate(result%filtered_covariance(2, 2, n))
    allocate(result%smoothed_state(2, n))
    allocate(result%smoothed_covariance(2, 2, n))
    result%filtered_state = filtered_mean
    result%filtered_covariance = filtered_covariance
    result%smoothed_state = smoothed_mean
    result%smoothed_covariance = smoothed_covariance
    result%transition_covariance = transition_covariance
    result%observation_variance = observation_variance
    result%initial_covariance = initial_covariance
  end subroutine dynamic_beta

  subroutine filter_and_smooth(y, x, initial_mean, initial_covariance, &
      transition_covariance, observation_variance, predicted_mean, &
      predicted_covariance, filtered_mean, filtered_covariance, &
      smoothed_mean, smoothed_covariance, smoother_gain)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: initial_mean(2)
    real(dp), intent(in) :: initial_covariance(2, 2)
    real(dp), intent(in) :: transition_covariance(2, 2)
    real(dp), intent(in) :: observation_variance
    real(dp), intent(out) :: predicted_mean(:, :)
    real(dp), intent(out) :: predicted_covariance(:, :, :)
    real(dp), intent(out) :: filtered_mean(:, :)
    real(dp), intent(out) :: filtered_covariance(:, :, :)
    real(dp), intent(out) :: smoothed_mean(:, :)
    real(dp), intent(out) :: smoothed_covariance(:, :, :)
    real(dp), intent(out) :: smoother_gain(:, :, :)
    real(dp) :: h(2)
    real(dp) :: gain(2)
    real(dp) :: innovation
    real(dp) :: innovation_variance
    real(dp) :: inv_pred(2, 2)
    real(dp) :: p_h(2)
    integer :: n
    integer :: t

    n = size(y)

    do t = 1, n
      if (t == 1) then
        predicted_mean(:, t) = initial_mean
        predicted_covariance(:, :, t) = initial_covariance
      else
        predicted_mean(:, t) = filtered_mean(:, t - 1)
        predicted_covariance(:, :, t) = &
          filtered_covariance(:, :, t - 1) + transition_covariance
      end if

      h = [x(t), 1.0_dp]
      p_h = matmul(predicted_covariance(:, :, t), h)
      innovation_variance = dot_product(h, p_h) + observation_variance
      innovation_variance = max(innovation_variance, 1.0e-12_dp)
      gain = p_h / innovation_variance
      innovation = y(t) - dot_product(h, predicted_mean(:, t))
      filtered_mean(:, t) = predicted_mean(:, t) + gain * innovation
      filtered_covariance(:, :, t) = predicted_covariance(:, :, t) - &
        outer_product(gain, p_h)
      call symmetrize_and_floor(filtered_covariance(:, :, t))
    end do

    smoothed_mean(:, n) = filtered_mean(:, n)
    smoothed_covariance(:, :, n) = filtered_covariance(:, :, n)

    do t = n - 1, 1, -1
      call inverse_2x2(predicted_covariance(:, :, t + 1), inv_pred)
      smoother_gain(:, :, t) = matmul(filtered_covariance(:, :, t), inv_pred)
      smoothed_mean(:, t) = filtered_mean(:, t) + &
        matmul(smoother_gain(:, :, t), &
        smoothed_mean(:, t + 1) - predicted_mean(:, t + 1))
      smoothed_covariance(:, :, t) = filtered_covariance(:, :, t) + &
        matmul(smoother_gain(:, :, t), matmul(&
        smoothed_covariance(:, :, t + 1) - predicted_covariance(:, :, t + 1), &
        transpose(smoother_gain(:, :, t))))
      call symmetrize_and_floor(smoothed_covariance(:, :, t))
    end do
  end subroutine filter_and_smooth

  subroutine em_update(y, x, initial_mean, initial_covariance, &
      transition_covariance, observation_variance, smoothed_mean, &
      smoothed_covariance, smoother_gain)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: initial_mean(2)
    real(dp), intent(inout) :: initial_covariance(2, 2)
    real(dp), intent(inout) :: transition_covariance(2, 2)
    real(dp), intent(inout) :: observation_variance
    real(dp), intent(in) :: smoothed_mean(:, :)
    real(dp), intent(in) :: smoothed_covariance(:, :, :)
    real(dp), intent(in) :: smoother_gain(:, :, :)
    real(dp) :: cross_covariance(2, 2)
    real(dp) :: difference(2)
    real(dp) :: h(2)
    real(dp) :: q_sum(2, 2)
    real(dp) :: residual
    real(dp) :: r_sum
    integer :: n
    integer :: t

    n = size(y)
    q_sum = 0.0_dp

    do t = 2, n
      difference = smoothed_mean(:, t) - smoothed_mean(:, t - 1)
      cross_covariance = matmul(smoothed_covariance(:, :, t), &
        transpose(smoother_gain(:, :, t - 1)))
      q_sum = q_sum + outer_product(difference, difference) + &
        smoothed_covariance(:, :, t) + smoothed_covariance(:, :, t - 1) - &
        cross_covariance - transpose(cross_covariance)
    end do
    transition_covariance = q_sum / real(n - 1, dp)
    call symmetrize_and_floor(transition_covariance)

    r_sum = 0.0_dp
    do t = 1, n
      h = [x(t), 1.0_dp]
      residual = y(t) - dot_product(h, smoothed_mean(:, t))
      r_sum = r_sum + residual**2 + &
        dot_product(h, matmul(smoothed_covariance(:, :, t), h))
    end do
    observation_variance = max(r_sum / real(n, dp), 1.0e-10_dp)

    difference = smoothed_mean(:, 1) - initial_mean
    initial_covariance = smoothed_covariance(:, :, 1) + &
      outer_product(difference, difference)
    call symmetrize_and_floor(initial_covariance)
  end subroutine em_update

  pure function outer_product(x, y) result(matrix)
    real(dp), intent(in) :: x(2)
    real(dp), intent(in) :: y(2)
    real(dp) :: matrix(2, 2)
    integer :: i
    integer :: j

    do j = 1, 2
      do i = 1, 2
        matrix(i, j) = x(i) * y(j)
      end do
    end do
  end function outer_product

  subroutine inverse_2x2(matrix, inverse)
    real(dp), intent(in) :: matrix(2, 2)
    real(dp), intent(out) :: inverse(2, 2)
    real(dp) :: determinant
    real(dp) :: work(2, 2)

    work = matrix
    determinant = work(1, 1) * work(2, 2) - work(1, 2) * work(2, 1)
    if (abs(determinant) < 1.0e-14_dp) then
      work(1, 1) = work(1, 1) + 1.0e-8_dp
      work(2, 2) = work(2, 2) + 1.0e-8_dp
      determinant = work(1, 1) * work(2, 2) - work(1, 2) * work(2, 1)
    end if

    inverse(1, 1) = work(2, 2) / determinant
    inverse(1, 2) = -work(1, 2) / determinant
    inverse(2, 1) = -work(2, 1) / determinant
    inverse(2, 2) = work(1, 1) / determinant
  end subroutine inverse_2x2

  subroutine symmetrize_and_floor(matrix)
    real(dp), intent(inout) :: matrix(2, 2)
    real(dp) :: off_diagonal

    off_diagonal = 0.5_dp * (matrix(1, 2) + matrix(2, 1))
    matrix(1, 2) = off_diagonal
    matrix(2, 1) = off_diagonal
    matrix(1, 1) = max(matrix(1, 1), 1.0e-12_dp)
    matrix(2, 2) = max(matrix(2, 2), 1.0e-12_dp)
  end subroutine symmetrize_and_floor

end module trading_dynamic_beta
