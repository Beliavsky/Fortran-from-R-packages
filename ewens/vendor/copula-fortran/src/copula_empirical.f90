! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_empirical
  use copula_kinds, only : dp, i8
  use copula_types, only : test_result
  use copula_random, only : seed_random, random_uniform
  implicit none
  private
  public :: pseudo_observations, empirical_copula, empirical_copula_values
  public :: sample_kendall_tau, sample_spearman_rho
  public :: independence_test, exchangeability_test, radial_symmetry_test
contains
  subroutine pseudo_observations(x, u, ok)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: u(:,:)
    logical, intent(out) :: ok
    integer :: n, d, i, j, k, less_count, equal_count
    n = size(x,1)
    d = size(x,2)
    allocate(u(n,d))
    ok = n >= 1 .and. d >= 1
    if (.not. ok) return
    do j = 1, d
      do i = 1, n
        less_count = 0
        equal_count = 0
        do k = 1, n
          if (x(k,j) < x(i,j)) then
            less_count = less_count+1
          else if (abs(x(k,j)-x(i,j)) <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i,j)))) then
            equal_count = equal_count+1
          end if
        end do
        u(i,j) = (real(less_count,dp)+0.5_dp*real(equal_count+1,dp))/real(n+1,dp)
      end do
    end do
  end subroutine pseudo_observations

  real(dp) function empirical_copula(u, point) result(value)
    real(dp), intent(in) :: u(:,:), point(:)
    integer :: i, n, count
    n = size(u,1)
    if (size(u,2) /= size(point) .or. n == 0) then
      value = 0.0_dp
      return
    end if
    count = 0
    do i = 1, n
      if (all(u(i,:) <= point)) count = count+1
    end do
    value = real(count,dp)/real(n,dp)
  end function empirical_copula

  subroutine empirical_copula_values(u, points, values, ok)
    real(dp), intent(in) :: u(:,:), points(:,:)
    real(dp), allocatable, intent(out) :: values(:)
    logical, intent(out) :: ok
    integer :: i
    ok = size(u,2) == size(points,2)
    allocate(values(size(points,1)))
    values = 0.0_dp
    if (.not. ok) return
    do i = 1, size(points,1)
      values(i) = empirical_copula(u,points(i,:))
    end do
  end subroutine empirical_copula_values

  real(dp) function sample_kendall_tau(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    integer :: i, j, n, concordant, discordant, ties_x, ties_y
    real(dp) :: denominator
    n = size(x)
    if (size(y) /= n .or. n < 2) then
      value = 0.0_dp
      return
    end if
    concordant = 0
    discordant = 0
    ties_x = 0
    ties_y = 0
    do i = 1, n-1
      do j = i+1, n
        if (abs(x(i)-x(j)) <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)))) ties_x = ties_x+1
        if (abs(y(i)-y(j)) <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(y(i)))) ties_y = ties_y+1
        if ((x(i)-x(j))*(y(i)-y(j)) > 0.0_dp) concordant = concordant+1
        if ((x(i)-x(j))*(y(i)-y(j)) < 0.0_dp) discordant = discordant+1
      end do
    end do
    denominator = sqrt(real(n*(n-1)/2-ties_x,dp)*real(n*(n-1)/2-ties_y,dp))
    if (denominator <= 0.0_dp) then
      value = 0.0_dp
    else
      value = real(concordant-discordant,dp)/denominator
    end if
  end function sample_kendall_tau

  real(dp) function sample_spearman_rho(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), allocatable :: data(:,:), ranks(:,:)
    logical :: ok
    integer :: n
    n = size(x)
    if (size(y) /= n .or. n < 2) then
      value = 0.0_dp
      return
    end if
    allocate(data(n,2))
    data(:,1) = x
    data(:,2) = y
    call pseudo_observations(data,ranks,ok)
    if (.not. ok) then
      value = 0.0_dp
    else
      value = correlation(ranks(:,1),ranks(:,2))
    end if
  end function sample_spearman_rho

  function independence_test(u, replicates, seed) result(result)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in), optional :: replicates
    integer(i8), intent(in), optional :: seed
    type(test_result) :: result
    real(dp), allocatable :: permuted(:,:)
    integer :: b, reps, exceed, j
    result%ok = size(u,1) >= 5 .and. size(u,2) >= 2
    if (.not. result%ok) then
      result%message = 'at least five multivariate observations are required'
      return
    end if
    reps = 199
    if (present(replicates)) reps = max(1,replicates)
    if (present(seed)) call seed_random(seed)
    result%statistic = independence_statistic(u)
    allocate(permuted(size(u,1),size(u,2)))
    exceed = 0
    do b = 1, reps
      permuted(:,1) = u(:,1)
      do j = 2, size(u,2)
        call permute_column(u(:,j),permuted(:,j))
      end do
      if (independence_statistic(permuted) >= result%statistic) exceed = exceed+1
    end do
    result%replicates = reps
    result%p_value = real(exceed+1,dp)/real(reps+1,dp)
  end function independence_test

  function exchangeability_test(u, replicates, seed) result(result)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in), optional :: replicates
    integer(i8), intent(in), optional :: seed
    type(test_result) :: result
    real(dp), allocatable :: work(:,:)
    integer :: reps, b, i, exceed
    result%ok = size(u,1) >= 5 .and. size(u,2) == 2
    if (.not. result%ok) then
      result%message = 'a bivariate sample with at least five rows is required'
      return
    end if
    reps = 199
    if (present(replicates)) reps = max(1,replicates)
    if (present(seed)) call seed_random(seed)
    result%statistic = exchangeability_statistic(u)
    allocate(work(size(u,1),2))
    exceed = 0
    do b = 1, reps
      work = u
      do i = 1, size(u,1)
        if (random_uniform() < 0.5_dp) work(i,:) = [u(i,2),u(i,1)]
      end do
      if (exchangeability_statistic(work) >= result%statistic) exceed = exceed+1
    end do
    result%replicates = reps
    result%p_value = real(exceed+1,dp)/real(reps+1,dp)
  end function exchangeability_test

  function radial_symmetry_test(u, replicates, seed) result(result)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in), optional :: replicates
    integer(i8), intent(in), optional :: seed
    type(test_result) :: result
    real(dp), allocatable :: work(:,:)
    integer :: reps, b, i, exceed
    result%ok = size(u,1) >= 5 .and. size(u,2) == 2
    if (.not. result%ok) then
      result%message = 'a bivariate sample with at least five rows is required'
      return
    end if
    reps = 199
    if (present(replicates)) reps = max(1,replicates)
    if (present(seed)) call seed_random(seed)
    result%statistic = radial_statistic(u)
    allocate(work(size(u,1),2))
    exceed = 0
    do b = 1, reps
      work = u
      do i = 1, size(u,1)
        if (random_uniform() < 0.5_dp) work(i,:) = 1.0_dp-u(i,:)
      end do
      if (radial_statistic(work) >= result%statistic) exceed = exceed+1
    end do
    result%replicates = reps
    result%p_value = real(exceed+1,dp)/real(reps+1,dp)
  end function radial_symmetry_test

  real(dp) function independence_statistic(u) result(value)
    real(dp), intent(in) :: u(:,:)
    integer :: i, n
    real(dp) :: difference
    n = size(u,1)
    value = 0.0_dp
    do i = 1, n
      difference = empirical_copula(u,u(i,:))-product(u(i,:))
      value = value+difference*difference
    end do
    value = value
  end function independence_statistic

  real(dp) function exchangeability_statistic(u) result(value)
    real(dp), intent(in) :: u(:,:)
    integer :: i, n
    real(dp) :: difference
    n = size(u,1)
    value = 0.0_dp
    do i = 1, n
      difference = empirical_copula(u,u(i,:))-empirical_copula(u,[u(i,2),u(i,1)])
      value = value+difference*difference
    end do
  end function exchangeability_statistic

  real(dp) function radial_statistic(u) result(value)
    real(dp), intent(in) :: u(:,:)
    real(dp), allocatable :: reflected(:,:)
    integer :: i, n
    real(dp) :: difference
    n = size(u,1)
    allocate(reflected(n,2))
    reflected = 1.0_dp-u
    value = 0.0_dp
    do i = 1, n
      difference = empirical_copula(u,u(i,:))-empirical_copula(reflected,u(i,:))
      value = value+difference*difference
    end do
  end function radial_statistic

  subroutine permute_column(input, output)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    integer :: i, j, n
    real(dp) :: temporary
    n = size(input)
    output = input
    do i = n, 2, -1
      j = 1+int(random_uniform()*real(i,dp))
      j = min(i,max(1,j))
      temporary = output(i)
      output(i) = output(j)
      output(j) = temporary
    end do
  end subroutine permute_column

  real(dp) function correlation(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: mx, my, sx, sy
    mx = sum(x)/real(size(x),dp)
    my = sum(y)/real(size(y),dp)
    sx = sqrt(sum((x-mx)**2))
    sy = sqrt(sum((y-my)**2))
    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = sum((x-mx)*(y-my))/(sx*sy)
    end if
  end function correlation
end module copula_empirical
