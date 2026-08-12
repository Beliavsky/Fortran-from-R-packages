module flsss_integerize
  use flsss_kinds, only : dp, i8
  use flsss_types, only : integerized_problem
  implicit none
  private
  public :: integerize_problem

contains

  function integerize_problem(len, v, target, me, precision_level) result(r)
    integer, intent(in) :: len
    real(dp), intent(in) :: v(:,:), target(:), me(:)
    integer, intent(in), optional :: precision_level(:)
    type(integerized_problem) :: r
    integer :: n, d, j, p, ratio, ratio_limit
    real(dp), allocatable :: vn(:)
    real(dp) :: the_max, least_diff

    n = size(v,1); d = size(v,2)
    if (len < 0 .or. len > n) error stop "integerize_problem: invalid subset size"
    ratio_limit = shiftr(huge(ratio), 1)
    if (size(target) /= d .or. size(me) /= d) error stop "integerize_problem: dimension mismatch"
    if (any(me <= 0.0_dp)) error stop "integerize_problem: ME must be positive"
    allocate(r%v(n,d), r%target(d), r%me(d), r%ratio(d), vn(n))

    do j = 1, d
      vn = v(:,j) / me(j)
      p = 0
      if (present(precision_level)) p = precision_level(j)
      if (p == 0) then
        ratio = 1
        the_max = max(0.0_dp, maxval(vn))
        do while (the_max * real(ratio,dp) < real(n*8,dp) .and. ratio < ratio_limit)
          ratio = ratio * 2
        end do
      else if (p == -1) then
        call minimum_positive_difference(vn, least_diff)
        ratio = 1
        if (least_diff > 0.0_dp) then
          do while (nint(least_diff * real(ratio,dp), kind=i8) < 1_i8 .and. ratio < ratio_limit)
            ratio = ratio * 2
          end do
        end if
      else
        ratio = 1
        the_max = max(0.0_dp, maxval(vn))
        do while (the_max * real(ratio,dp) < real(max(1,p),dp) .and. ratio < ratio_limit)
          ratio = ratio * 2
        end do
      end if
      r%ratio(j) = ratio
      r%v(:,j) = nint(vn * real(ratio,dp), kind=i8)
      r%target(j) = nint(target(j) / me(j) * real(ratio,dp), kind=i8)
      r%me(j) = int(ratio, i8)
    end do
    r%compressed_dim = d
  end function integerize_problem

  subroutine minimum_positive_difference(x, diffmin)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: diffmin
    real(dp), allocatable :: y(:)
    real(dp) :: t
    integer :: i, j
    allocate(y(size(x))); y = x
    do i = 2, size(y)
      t = y(i); j = i - 1
      do while (j >= 1)
        if (y(j) <= t) exit
        y(j+1) = y(j); j = j - 1
      end do
      y(j+1) = t
    end do
    diffmin = huge(1.0_dp)
    do i = 2, size(y)
      t = y(i) - y(i-1)
      if (t > 1.0e-10_dp .and. t < diffmin) diffmin = t
    end do
    if (diffmin > 0.5_dp * huge(1.0_dp)) diffmin = 0.0_dp
  end subroutine minimum_positive_difference

end module flsss_integerize
