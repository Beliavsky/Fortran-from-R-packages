module rmalschains_operators
  use rmalschains_kinds, only : dp
  use rmalschains_rng, only : rng_state, rng_uniform, rng_int
  implicit none
  private
  public :: blx_alpha, bga_mutate, nam_select, nearest_distance, nearest_coordinate_distances
contains
  subroutine blx_alpha(mom, dad, lower, upper, alpha, rng, child)
    real(dp), intent(in) :: mom(:), dad(:), lower(:), upper(:), alpha
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: child(:)
    integer :: i
    real(dp) :: x1, x2, ext, a, b
    if (size(child) /= size(mom) .or. size(dad) /= size(mom)) error stop "blx_alpha: size mismatch"
    do i = 1, size(mom)
      x1 = min(mom(i), dad(i))
      x2 = max(mom(i), dad(i))
      ext = alpha * (x2 - x1)
      a = max(lower(i), x1 - ext)
      b = min(upper(i), x2 + ext)
      child(i) = a + rng_uniform(rng) * (b - a)
    end do
  end subroutine blx_alpha

  subroutine bga_mutate(x, lower, upper, rng, changed)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    type(rng_state), intent(inout) :: rng
    logical, intent(out), optional :: changed
    integer :: pos, i
    real(dp) :: diff, s, val
    logical :: did
    did = rng_uniform(rng) <= 0.125_dp
    if (did) then
      pos = rng_int(rng, 1, size(x))
      s = 0.0_dp
      diff = 1.0_dp
      do i = 1, 16
        if (rng_uniform(rng) < 1.0_dp / 16.0_dp) s = s + diff
        diff = diff / 2.0_dp
      end do
      if (s > 0.0_dp) then
        val = x(pos)
        if (rng_uniform(rng) < 0.5_dp) then
          val = val + 0.1_dp * (upper(pos) - lower(pos)) * s
        else
          val = val - 0.1_dp * (upper(pos) - lower(pos)) * s
        end if
        x(pos) = min(max(val, lower(pos)), upper(pos))
      end if
    end if
    if (present(changed)) changed = did
  end subroutine bga_mutate

  subroutine nam_select(pop, rng, mom, dad)
    real(dp), intent(in) :: pop(:, :)
    type(rng_state), intent(inout) :: rng
    integer, intent(out) :: mom, dad
    integer :: n, k, cand, tries
    logical, allocatable :: used(:)
    real(dp) :: dist, best_dist
    n = size(pop, 2)
    if (n < 4) error stop "nam_select: population must contain at least four individuals"
    allocate(used(n))
    used = .false.
    mom = rng_int(rng, 1, n)
    used(mom) = .true.
    best_dist = -1.0_dp
    dad = 0
    do k = 1, 3
      tries = 0
      do
        cand = rng_int(rng, 1, n)
        tries = tries + 1
        if (.not. used(cand)) exit
        if (tries > 10 * n) error stop "nam_select: sampling failure"
      end do
      used(cand) = .true.
      dist = sqrt(sum((pop(:, mom) - pop(:, cand))**2)) / real(size(pop, 1), dp)
      if (dist > best_dist) then
        best_dist = dist
        dad = cand
      end if
    end do
  end subroutine nam_select

  function nearest_distance(x, pop, skip) result(dist)
    real(dp), intent(in) :: x(:), pop(:, :)
    integer, intent(in), optional :: skip
    real(dp) :: dist
    integer :: j
    real(dp) :: d
    dist = huge(1.0_dp)
    do j = 1, size(pop, 2)
      if (present(skip)) then
        if (j == skip) cycle
      end if
      d = sqrt(sum((x - pop(:, j))**2)) / real(size(x), dp)
      if (d > 0.0_dp .and. d < dist) dist = d
    end do
    if (dist >= huge(1.0_dp) / 2.0_dp) dist = 0.0_dp
  end function nearest_distance

  subroutine nearest_coordinate_distances(x, pop, dist, skip)
    real(dp), intent(in) :: x(:), pop(:, :)
    real(dp), intent(out) :: dist(:)
    integer, intent(in), optional :: skip
    integer :: i, j
    real(dp) :: d
    if (size(dist) /= size(x)) error stop "nearest_coordinate_distances: size mismatch"
    dist = huge(1.0_dp)
    do j = 1, size(pop, 2)
      if (present(skip)) then
        if (j == skip) cycle
      end if
      do i = 1, size(x)
        d = abs(x(i) - pop(i, j))
        if (d > 0.0_dp .and. d < dist(i)) dist(i) = d
      end do
    end do
    do i = 1, size(x)
      if (dist(i) >= huge(1.0_dp) / 2.0_dp) dist(i) = 0.0_dp
    end do
  end subroutine nearest_coordinate_distances
end module rmalschains_operators
