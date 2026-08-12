module flsss_util
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions, subset_solution
  implicit none
  private
  public :: timer_type, append_solution, argsort_real, top_k_sum

  type :: timer_type
    integer(i8) :: start_count = 0_i8
    integer(i8) :: rate = 1_i8
    real(dp) :: limit = huge(1.0_dp)
  contains
    procedure :: start => timer_start
    procedure :: expired => timer_expired
    procedure :: elapsed => timer_elapsed
  end type timer_type

contains

  subroutine timer_start(this, limit)
    class(timer_type), intent(inout) :: this
    real(dp), intent(in), optional :: limit
    integer(i8) :: mx
    call system_clock(this%start_count, this%rate, mx)
    if (present(limit)) this%limit = max(0.0_dp, limit)
  end subroutine timer_start

  logical function timer_expired(this) result(tf)
    class(timer_type), intent(in) :: this
    tf = this%elapsed() >= this%limit
  end function timer_expired

  real(dp) function timer_elapsed(this) result(t)
    class(timer_type), intent(in) :: this
    integer(i8) :: now
    call system_clock(now)
    t = real(now - this%start_count, dp) / real(max(1_i8, this%rate), dp)
  end function timer_elapsed

  subroutine append_solution(r, idx)
    type(subset_solutions), intent(inout) :: r
    integer, intent(in) :: idx(:)
    type(subset_solution), allocatable :: tmp(:)
    integer :: n
    if (.not. allocated(r%sol)) then
      allocate(r%sol(1))
      allocate(r%sol(1)%idx(size(idx)))
      r%sol(1)%idx = idx
      return
    end if
    n = size(r%sol)
    allocate(tmp(n + 1))
    tmp(1:n) = r%sol
    allocate(tmp(n + 1)%idx(size(idx)))
    tmp(n + 1)%idx = idx
    call move_alloc(tmp, r%sol)
  end subroutine append_solution

  function argsort_real(x, descending) result(ord)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: descending
    integer :: ord(size(x)), i, j, key
    logical :: desc
    desc = .false.
    if (present(descending)) desc = descending
    do i = 1, size(x)
      ord(i) = i
    end do
    do i = 2, size(x)
      key = ord(i)
      j = i - 1
      if (desc) then
        do while (j >= 1)
          if (x(ord(j)) >= x(key)) exit
          ord(j + 1) = ord(j)
          j = j - 1
        end do
      else
        do while (j >= 1)
          if (x(ord(j)) <= x(key)) exit
          ord(j + 1) = ord(j)
          j = j - 1
        end do
      end if
      ord(j + 1) = key
    end do
  end function argsort_real

  real(dp) function top_k_sum(x, k) result(s)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: k
    integer :: o(size(x)), kk
    kk = min(max(k, 0), size(x))
    if (kk == 0) then
      s = 0.0_dp
      return
    end if
    o = argsort_real(x, .true.)
    s = sum(x(o(1:kk)))
  end function top_k_sum

end module flsss_util
