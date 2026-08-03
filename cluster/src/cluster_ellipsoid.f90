! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_ellipsoid
  use fastcluster_kinds, only: dp
  use cluster_types, only: ellipsoid_result, cluster_success, cluster_invalid_argument, &
    cluster_numerical_failure, cluster_not_converged
  use cluster_linalg, only: inverse_matrix, determinant_spd, symmetric_eigen
  implicit none
  private

  public :: ellipsoidhull
  public :: predict_ellipsoid
  public :: volume_ellipsoid
  public :: ellipsoid_points

contains

  subroutine ellipsoidhull(x, result, tolerance, max_iter)
    real(dp), intent(in) :: x(:, :)
    type(ellipsoid_result), intent(out) :: result
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iter

    real(dp), allocatable :: q(:, :), weighted(:, :), inverse(:, :), u(:), unew(:), centered(:, :)
    real(dp) :: max_m, max_quad, step, tol, determinant
    integer :: d, i, imax, iter, n, niter, status

    n = size(x, 1)
    d = size(x, 2)
    tol = 1.0e-7_dp
    if (present(tolerance)) tol = tolerance
    niter = 5000
    if (present(max_iter)) niter = max_iter
    if (n <= d .or. d < 1 .or. tol <= 0.0_dp .or. niter < 1) then
      call fail_ellipsoid(result, cluster_invalid_argument, &
        'ellipsoidhull requires n > dimension and positive controls')
      return
    end if
    allocate(q(d+1, n), weighted(d+1, d+1), u(n), unew(n), centered(n, d))
    q(1:d, :) = transpose(x)
    q(d+1, :) = 1.0_dp
    u = 1.0_dp / real(n, dp)

    do iter = 1, niter
      weighted = 0.0_dp
      do i = 1, n
        weighted = weighted + u(i) * outer_product(q(:, i), q(:, i))
      end do
      call inverse_matrix(weighted, inverse, status)
      if (status /= 0) then
        call fail_ellipsoid(result, cluster_numerical_failure, &
          'ellipsoid iteration encountered a singular matrix')
        return
      end if
      max_m = -huge(1.0_dp)
      imax = 1
      do i = 1, n
        step = dot_product(q(:, i), matmul(inverse, q(:, i)))
        if (step > max_m) then
          max_m = step
          imax = i
        end if
      end do
      if (max_m <= real(d+1, dp) * (1.0_dp + tol)) exit
      step = (max_m - real(d+1, dp)) / &
        (real(d+1, dp) * (max_m - 1.0_dp))
      step = min(1.0_dp, max(0.0_dp, step))
      unew = (1.0_dp - step) * u
      unew(imax) = unew(imax) + step
      if (maxval(abs(unew - u)) <= tol) then
        u = unew
        exit
      end if
      u = unew
    end do

    allocate(result%center(d), result%covariance(d, d), result%shape(d, d))
    result%center = matmul(transpose(x), u)
    do i = 1, n
      centered(i, :) = x(i, :) - result%center
    end do
    result%covariance = 0.0_dp
    do i = 1, n
      result%covariance = result%covariance + u(i) * &
        outer_product(centered(i, :), centered(i, :))
    end do
    result%covariance = real(d, dp) * result%covariance
    call inverse_matrix(result%covariance, result%shape, status)
    if (status /= 0) then
      call fail_ellipsoid(result, cluster_numerical_failure, &
        'final ellipsoid covariance is singular')
      return
    end if
    max_quad = 0.0_dp
    do i = 1, n
      max_quad = max(max_quad, dot_product(centered(i, :), &
        matmul(result%shape, centered(i, :))))
    end do
    if (max_quad > 1.0_dp) then
      result%covariance = max_quad * result%covariance
      result%shape = result%shape / max_quad
    end if
    call determinant_spd(result%covariance, determinant, status)
    if (status /= 0) then
      call fail_ellipsoid(result, cluster_numerical_failure, &
        'final ellipsoid covariance is not positive definite')
      return
    end if
    result%radius2 = 1.0_dp
    result%volume = unit_ball_volume(d) * sqrt(determinant)
    result%iterations = min(iter, niter)
    result%status = cluster_success
    if (iter > niter) then
      result%message = 'ok; iteration cap reached and enclosure rescaled'
    else
      result%message = 'ok'
    end if
  end subroutine ellipsoidhull

  subroutine predict_ellipsoid(result, points, squared_distance, inside, status)
    type(ellipsoid_result), intent(in) :: result
    real(dp), intent(in) :: points(:, :)
    real(dp), allocatable, intent(out) :: squared_distance(:)
    logical, allocatable, intent(out), optional :: inside(:)
    integer, intent(out), optional :: status

    real(dp), allocatable :: delta(:)
    integer :: i, n

    n = size(points, 1)
    if (.not. result%ok() .or. size(points, 2) /= size(result%center)) then
      allocate(squared_distance(0))
      if (present(inside)) allocate(inside(0))
      if (present(status)) status = cluster_invalid_argument
      return
    end if
    allocate(squared_distance(n), delta(size(result%center)))
    if (present(inside)) allocate(inside(n))
    do i = 1, n
      delta = points(i, :) - result%center
      squared_distance(i) = dot_product(delta, matmul(result%shape, delta))
      if (present(inside)) inside(i) = squared_distance(i) <= result%radius2 * &
        (1.0_dp + 100.0_dp*epsilon(1.0_dp))
    end do
    if (present(status)) status = cluster_success
  end subroutine predict_ellipsoid

  real(dp) function volume_ellipsoid(result) result(value)
    type(ellipsoid_result), intent(in) :: result
    if (result%ok()) then
      value = result%volume
    else
      value = 0.0_dp
    end if
  end function volume_ellipsoid

  subroutine ellipsoid_points(result, n_points, points, status)
    type(ellipsoid_result), intent(in) :: result
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: points(:, :)
    integer, intent(out), optional :: status

    real(dp), allocatable :: values(:), vectors(:, :), transform(:, :), direction(:)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: angle, phi, z
    integer :: d, i, st

    if (.not. result%ok() .or. n_points < 1) then
      allocate(points(0, 0))
      if (present(status)) status = cluster_invalid_argument
      return
    end if
    d = size(result%center)
    if (d /= 2 .and. d /= 3) then
      allocate(points(0, 0))
      if (present(status)) status = cluster_invalid_argument
      return
    end if
    call symmetric_eigen(result%covariance, values, vectors, st)
    if (st /= 0 .or. any(values <= 0.0_dp)) then
      allocate(points(0, 0))
      if (present(status)) status = cluster_numerical_failure
      return
    end if
    allocate(transform(d, d), direction(d), points(n_points, d))
    transform = matmul(vectors, diagonal_matrix(sqrt(values)))
    do i = 1, n_points
      if (d == 2) then
        angle = 2.0_dp*pi*real(i-1, dp)/real(n_points, dp)
        direction = [cos(angle), sin(angle)]
      else
        z = 1.0_dp - 2.0_dp*(real(i, dp)-0.5_dp)/real(n_points, dp)
        phi = pi*(3.0_dp-sqrt(5.0_dp))*real(i-1, dp)
        direction = [sqrt(max(0.0_dp, 1.0_dp-z*z))*cos(phi), &
          sqrt(max(0.0_dp, 1.0_dp-z*z))*sin(phi), z]
      end if
      points(i, :) = result%center + matmul(transform, direction)
    end do
    if (present(status)) status = cluster_success
  end subroutine ellipsoid_points

  pure function outer_product(a, b) result(output)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: output(size(a), size(b))
    integer :: i
    do i = 1, size(a)
      output(i, :) = a(i) * b
    end do
  end function outer_product

  pure function diagonal_matrix(values) result(matrix)
    real(dp), intent(in) :: values(:)
    real(dp) :: matrix(size(values), size(values))
    integer :: i
    matrix = 0.0_dp
    do i = 1, size(values)
      matrix(i, i) = values(i)
    end do
  end function diagonal_matrix

  real(dp) function unit_ball_volume(d) result(value)
    integer, intent(in) :: d
    real(dp), parameter :: pi = acos(-1.0_dp)
    value = pi**(0.5_dp*real(d, dp)) / gamma(0.5_dp*real(d, dp) + 1.0_dp)
  end function unit_ball_volume

  subroutine fail_ellipsoid(result, status, message)
    type(ellipsoid_result), intent(out) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_ellipsoid

end module cluster_ellipsoid
