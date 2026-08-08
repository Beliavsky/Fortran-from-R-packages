! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_operators
  use rgenoud_kinds, only : dp
  use rgenoud_random, only : randu, randi, coin_flip
  implicit none
  private
  integer, parameter :: max_unique_try = 1000
  public :: oper_uniform, oper_boundary, oper_nonuniform, oper_polytope
  public :: oper_simple, oper_whole_nonuniform, oper_heuristic, in_bounds
contains
  logical function in_bounds(x, lower, upper) result(ok)
    real(dp), intent(in) :: x(:), lower(:), upper(:)
    ok = all(x >= lower .and. x <= upper)
  end function in_bounds

  subroutine find_range(parent, lower, upper, comp, lo, hi)
    real(dp), intent(in) :: parent(:), lower(:), upper(:)
    integer, intent(in) :: comp
    real(dp), intent(out) :: lo, hi
    real(dp) :: a
    a = randu()
    lo = a * lower(comp) + (1.0_dp - a) * parent(comp)
    a = randu()
    hi = (1.0_dp - a) * parent(comp) + a * upper(comp)
  end subroutine find_range

  real(dp) function get_f(total_gen, gen, y, b) result(delta)
    integer, intent(in) :: total_gen, gen, b
    real(dp), intent(in) :: y
    real(dp) :: factor
    factor = (1.0_dp - real(gen, dp) / real(max(1, total_gen), dp))**real(b, dp)
    factor = factor * randu()
    factor = max(factor, 1.0e-5_dp)
    delta = y * factor
  end function get_f

  subroutine oper_uniform(parent, lower, upper, is_integer)
    real(dp), intent(inout) :: parent(:)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(in) :: is_integer
    integer :: comp, it, ilo, ihi, iv
    real(dp) :: lo, hi, value
    do it = 1, max_unique_try
      comp = randi(1, size(parent))
      call find_range(parent, lower, upper, comp, lo, hi)
      if (is_integer) then
        ilo = max(ceiling(lo), ceiling(lower(comp)))
        ihi = min(floor(hi), floor(upper(comp)))
        if (ihi < ilo) cycle
        iv = randi(ilo, ihi)
        value = real(iv, dp)
      else
        value = randu(lo, hi)
      end if
      if (abs(value - parent(comp)) > 0.0_dp .or. it == max_unique_try) then
        parent(comp) = value
        return
      end if
    end do
  end subroutine oper_uniform

  subroutine oper_boundary(parent, lower, upper, is_integer)
    real(dp), intent(inout) :: parent(:)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(in) :: is_integer
    integer :: comp, it, ilo, ihi
    real(dp) :: lo, hi, value
    do it = 1, max_unique_try
      comp = randi(1, size(parent))
      call find_range(parent, lower, upper, comp, lo, hi)
      if (is_integer) then
        ilo = max(ceiling(lo), ceiling(lower(comp)))
        ihi = min(floor(hi), floor(upper(comp)))
        value = real(merge(ihi, ilo, coin_flip()), dp)
      else
        value = merge(hi, lo, coin_flip())
      end if
      if (abs(value - parent(comp)) > 0.0_dp .or. it == max_unique_try) then
        parent(comp) = value
        return
      end if
    end do
  end subroutine oper_boundary

  subroutine oper_nonuniform(parent, lower, upper, total_gen, gen, b, is_integer)
    real(dp), intent(inout) :: parent(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: total_gen, gen, b
    logical, intent(in) :: is_integer
    integer :: comp, it
    real(dp) :: lo, hi, value
    do it = 1, max_unique_try
      comp = randi(1, size(parent))
      call find_range(parent, lower, upper, comp, lo, hi)
      if (coin_flip()) then
        value = parent(comp) + get_f(total_gen, gen, hi - parent(comp), b)
      else
        value = parent(comp) - get_f(total_gen, gen, parent(comp) - lo, b)
      end if
      if (is_integer) value = real(int(value), dp)
      value = min(max(value, lower(comp)), upper(comp))
      if (abs(value - parent(comp)) > 0.0_dp .or. it == max_unique_try) then
        parent(comp) = value
        return
      end if
    end do
  end subroutine oper_nonuniform

  subroutine oper_polytope(parents, child, lower, upper, is_integer)
    real(dp), intent(in) :: parents(:, :)
    real(dp), intent(out) :: child(:)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(in) :: is_integer
    real(dp) :: w(size(parents, 2)), sw
    integer :: k
    do k = 1, size(w)
      do
        w(k) = randu()
        if (w(k) > 0.0_dp) exit
      end do
    end do
    sw = sum(w)
    w = w / sw
    child = matmul(parents, w)
    if (is_integer) child = real(int(child), dp)
    child = min(max(child, lower), upper)
  end subroutine oper_polytope

  subroutine oper_simple(p1, p2, lower, upper, is_integer)
    real(dp), intent(inout) :: p1(:), p2(:)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(in) :: is_integer
    real(dp) :: c1(size(p1)), c2(size(p2)), a
    integer :: cut, n, step, i, it
    step = 10
    do it = 1, max_unique_try
      cut = randi(1, size(p1))
      c1 = p1
      c2 = p2
      do n = 1, step
        a = real(n, dp) / real(step, dp)
        do i = cut + 1, size(p1)
          c1(i) = a * p1(i) + (1.0_dp - a) * p2(i)
          c2(i) = a * p2(i) + (1.0_dp - a) * p1(i)
        end do
        if (is_integer) then
          c1 = real(int(c1), dp)
          c2 = real(int(c2), dp)
        end if
        if (in_bounds(c1, lower, upper) .and. in_bounds(c2, lower, upper)) exit
      end do
      if (any(abs(c1 - p1) > 0.0_dp) .and. any(abs(c2 - p2) > 0.0_dp)) then
        p1 = c1
        p2 = c2
        return
      end if
    end do
  end subroutine oper_simple

  subroutine oper_whole_nonuniform(parent, lower, upper, total_gen, gen, b, is_integer)
    real(dp), intent(inout) :: parent(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: total_gen, gen, b
    logical, intent(in) :: is_integer
    real(dp) :: lo, hi, value
    integer :: i
    do i = 1, size(parent)
      call find_range(parent, lower, upper, i, lo, hi)
      if (coin_flip()) then
        value = parent(i) + get_f(total_gen, gen, hi - parent(i), b)
      else
        value = parent(i) - get_f(total_gen, gen, parent(i) - lo, b)
      end if
      if (is_integer) value = real(int(value), dp)
      parent(i) = min(max(value, lower(i)), upper(i))
    end do
  end subroutine oper_whole_nonuniform

  subroutine oper_heuristic(worse, better, lower, upper, is_integer)
    real(dp), intent(inout) :: worse(:)
    real(dp), intent(in) :: better(:), lower(:), upper(:)
    logical, intent(in) :: is_integer
    real(dp) :: child(size(worse)), a
    integer :: it
    do it = 1, max_unique_try
      a = randu()
      child = a * (better - worse) + better
      if (is_integer) child = real(int(child), dp)
      if (in_bounds(child, lower, upper) .and. any(abs(child - worse) > 0.0_dp)) then
        worse = child
        return
      end if
    end do
  end subroutine oper_heuristic
end module rgenoud_operators
