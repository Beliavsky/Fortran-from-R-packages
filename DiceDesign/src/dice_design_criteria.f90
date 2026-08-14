module dice_design_criteria
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  use dice_design_kinds, only : dp
  use dice_design_utils, only : pairwise_distances, rescale_unit_cube, mean_real, sample_sd
  implicit none
  private

  type, public :: mst_result
    integer, allocatable :: edges(:, :)
    real(dp), allocatable :: lengths(:)
    real(dp) :: mean_length = 0.0_dp
    real(dp) :: sd_length = 0.0_dp
  end type mst_result

  public :: mindist, phi_p, mesh_ratio, coverage, discrepancy_all, discrepancy_value
  public :: mst_criteria

contains

  function mindist(design) result(value)
    real(dp), intent(in) :: design(:, :)
    real(dp) :: value
    real(dp), allocatable :: x(:, :), dmat(:, :)
    logical :: changed
    integer :: n, i, j

    call rescale_unit_cube(design, x, changed)
    n = size(x, 1)
    if (n < 2) then
      value = ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if
    call pairwise_distances(x, dmat)
    value = huge(1.0_dp)
    do j = 2, n
      do i = 1, j - 1
        value = min(value, dmat(i, j))
      end do
    end do
  end function mindist

  function phi_p(design, p) result(value)
    real(dp), intent(in) :: design(:, :)
    real(dp), intent(in), optional :: p
    real(dp) :: value, pp, d
    integer :: n, i, j

    pp = 50.0_dp
    if (present(p)) pp = p
    n = size(design, 1)
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    value = 0.0_dp
    do j = 2, n
      do i = 1, j - 1
        d = sqrt(sum((design(i, :) - design(j, :))**2))
        if (d <= 0.0_dp) then
          value = ieee_value(0.0_dp, ieee_positive_inf)
          return
        end if
        value = value + d**(-pp)
      end do
    end do
    value = value**(1.0_dp / pp)
  end function phi_p

  function mesh_ratio(design) result(value)
    real(dp), intent(in) :: design(:, :)
    real(dp) :: value
    real(dp), allocatable :: x(:, :)
    real(dp) :: d2, nearest, min_nearest, max_nearest
    logical :: changed
    integer :: n, d, i, j

    call rescale_unit_cube(design, x, changed)
    n = size(x, 1)
    d = size(x, 2)
    if (n < d) error stop 'mesh_ratio: number of points is lower than dimension'
    if (n < 2) then
      value = 1.0_dp
      return
    end if
    min_nearest = huge(1.0_dp)
    max_nearest = 0.0_dp
    do i = 1, n
      nearest = huge(1.0_dp)
      do j = 1, n
        if (i == j) cycle
        d2 = sum((x(i, :) - x(j, :))**2)
        nearest = min(nearest, d2)
      end do
      min_nearest = min(min_nearest, nearest)
      max_nearest = max(max_nearest, nearest)
    end do
    if (min_nearest <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_positive_inf)
    else
      value = sqrt(max_nearest / min_nearest)
    end if
  end function mesh_ratio

  function coverage(design) result(value)
    real(dp), intent(in) :: design(:, :)
    real(dp) :: value
    real(dp), allocatable :: x(:, :), nearest(:)
    real(dp) :: d2, gamma_bar
    logical :: changed
    integer :: n, d, i, j

    call rescale_unit_cube(design, x, changed)
    n = size(x, 1)
    d = size(x, 2)
    if (n < d) error stop 'coverage: number of points is lower than dimension'
    if (n < 2) then
      value = 0.0_dp
      return
    end if
    allocate(nearest(n))
    do i = 1, n
      nearest(i) = huge(1.0_dp)
      do j = 1, n
        if (i == j) cycle
        d2 = sum((x(i, :) - x(j, :))**2)
        nearest(i) = min(nearest(i), sqrt(d2))
      end do
    end do
    gamma_bar = mean_real(nearest)
    if (gamma_bar <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_positive_inf)
    else
      value = sqrt(sum((nearest - gamma_bar)**2) / real(n, dp)) / gamma_bar
    end if
  end function coverage

  subroutine discrepancy_all(design, values)
    real(dp), intent(in) :: design(:, :)
    real(dp), intent(out) :: values(7)
    real(dp), allocatable :: x(:, :)
    logical :: changed
    real(dp) :: s1, s2, p, q, t, t1, t2, dx
    integer :: n, d, i, k, j

    call rescale_unit_cube(design, x, changed)
    n = size(x, 1)
    d = size(x, 2)
    if (n < d) error stop 'discrepancy_all: number of points is lower than dimension'

    ! Centered L2 discrepancy.
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      p = product(1.0_dp + 0.5_dp*abs(x(i, :) - 0.5_dp) - &
                  0.5_dp*abs(x(i, :) - 0.5_dp)**2)
      s1 = s1 + p
      do k = 1, n
        q = product(1.0_dp + 0.5_dp*abs(x(i, :) - 0.5_dp) + &
                    0.5_dp*abs(x(k, :) - 0.5_dp) - 0.5_dp*abs(x(i, :) - x(k, :)))
        s2 = s2 + q
      end do
    end do
    values(1) = sqrt(max(0.0_dp, (13.0_dp/12.0_dp)**d - 2.0_dp*s1/real(n,dp) + s2/real(n*n,dp)))

    ! L2 discrepancy.
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      s1 = s1 + product(x(i, :) * (1.0_dp - x(i, :)))
      do k = 1, n
        q = 1.0_dp
        do j = 1, d
          q = q * (1.0_dp - max(x(i,j), x(k,j))) * min(x(i,j), x(k,j))
        end do
        s2 = s2 + q
      end do
    end do
    values(2) = sqrt(max(0.0_dp, 12.0_dp**(-d) - &
      (2.0_dp**(1-d)/real(n,dp))*s1 + s2/real(n*n,dp)))

    ! Star L2 discrepancy.
    t = 0.0_dp
    do j = 1, n
      do i = 1, n
        if (i /= j) then
          q = product(1.0_dp - max(x(i, :), x(j, :))) / real(n*n, dp)
        else
          t1 = product(1.0_dp - x(i, :))
          t2 = product(1.0_dp - x(i, :)**2)
          q = t1/real(n*n,dp) - (2.0_dp**(1-d)/real(n,dp))*t2
        end if
        t = t + q
      end do
    end do
    values(3) = sqrt(max(0.0_dp, 3.0_dp**(-d) + t))

    ! Modified L2 discrepancy.
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      s1 = s1 + product(3.0_dp - x(i, :)**2)
      do k = 1, n
        q = 1.0_dp
        do j = 1, d
          q = q * (2.0_dp - max(x(i,j), x(k,j)))
        end do
        s2 = s2 + q
      end do
    end do
    values(4) = sqrt(max(0.0_dp, (4.0_dp/3.0_dp)**d - &
      (2.0_dp**(1-d)/real(n,dp))*s1 + s2/real(n*n,dp)))

    ! Symmetric L2 discrepancy.
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      s1 = s1 + product(1.0_dp + 2.0_dp*x(i, :) - 2.0_dp*x(i, :)**2)
      do k = 1, n
        s2 = s2 + product(1.0_dp - abs(x(i, :) - x(k, :)))
      end do
    end do
    values(5) = sqrt(max(0.0_dp, (4.0_dp/3.0_dp)**d - 2.0_dp*s1/real(n,dp) + &
      2.0_dp**d*s2/real(n*n,dp)))

    ! Wrap-around L2 discrepancy.
    s1 = 0.0_dp
    do i = 1, n
      do k = 1, n
        q = 1.0_dp
        do j = 1, d
          dx = abs(x(i,j)-x(k,j))
          q = q * (1.5_dp - dx*(1.0_dp-dx))
        end do
        s1 = s1 + q
      end do
    end do
    values(6) = sqrt(max(0.0_dp, -(4.0_dp/3.0_dp)**d + s1/real(n*n,dp)))

    ! Mixture L2 discrepancy.
    s1 = 0.0_dp
    s2 = 0.0_dp
    do i = 1, n
      s1 = s1 + product(5.0_dp/3.0_dp - 0.25_dp*abs(x(i,:)-0.5_dp) - &
        0.25_dp*abs(x(i,:)-0.5_dp)**2)
      do k = 1, n
        s2 = s2 + product(15.0_dp/8.0_dp - 0.25_dp*abs(x(i,:)-0.5_dp) - &
          0.25_dp*abs(x(k,:)-0.5_dp) - 0.75_dp*abs(x(i,:)-x(k,:)) + &
          0.5_dp*abs(x(i,:)-x(k,:))**2)
      end do
    end do
    values(7) = sqrt(max(0.0_dp, (19.0_dp/12.0_dp)**d - 2.0_dp*s1/real(n,dp) + &
      s2/real(n*n,dp)))
  end subroutine discrepancy_all

  function discrepancy_value(design, kind) result(value)
    real(dp), intent(in) :: design(:, :)
    character(len=*), intent(in) :: kind
    real(dp) :: value, vals(7)
    character(len=:), allocatable :: k

    call discrepancy_all(design, vals)
    k = trim(adjustl(kind))
    select case (k)
    case ('C2', 'c2')
      value = vals(1)
    case ('L2', 'l2')
      value = vals(2)
    case ('L2star', 'l2star', 'L2STAR')
      value = vals(3)
    case ('M2', 'm2')
      value = vals(4)
    case ('S2', 's2')
      value = vals(5)
    case ('W2', 'w2')
      value = vals(6)
    case ('Mix2', 'mix2', 'MIX2')
      value = vals(7)
    case default
      error stop 'discrepancy_value: unknown discrepancy type'
    end select
  end function discrepancy_value

  subroutine mst_criteria(design, result)
    real(dp), intent(in) :: design(:, :)
    type(mst_result), intent(out) :: result
    real(dp), allocatable :: dmat(:, :)
    logical, allocatable :: in_tree(:)
    real(dp) :: best
    integer :: n, e, i, j, bi, bj

    n = size(design, 1)
    if (n <= 1) then
      allocate(result%edges(2, 0), result%lengths(0))
      return
    end if
    call pairwise_distances(design, dmat)
    allocate(in_tree(n), result%edges(2, n-1), result%lengths(n-1))
    in_tree = .false.
    in_tree(1) = .true.
    do e = 1, n - 1
      best = huge(1.0_dp)
      bi = 0
      bj = 0
      do i = 1, n
        if (.not. in_tree(i)) cycle
        do j = 1, n
          if (in_tree(j)) cycle
          if (dmat(i,j) < best) then
            best = dmat(i,j)
            bi = i
            bj = j
          end if
        end do
      end do
      if (bj == 0) error stop 'mst_criteria: disconnected distance graph'
      result%edges(:,e) = [bi, bj]
      result%lengths(e) = best
      in_tree(bj) = .true.
    end do
    result%mean_length = mean_real(result%lengths)
    result%sd_length = sample_sd(result%lengths)
  end subroutine mst_criteria

end module dice_design_criteria
