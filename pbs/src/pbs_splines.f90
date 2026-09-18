module pbs_splines
   use pbs_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_next_after
   implicit none
   private

   public :: bspline_design
   public :: quantile_type7
   public :: sort_real

contains

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Values to sort in ascending order in place.
      integer :: i
      integer :: j
      real(dp) :: key

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real

   pure function quantile_type7(x, prob) result(q)
      real(dp), intent(in) :: x(:) !! Finite sample values whose type-7 quantile is required.
      real(dp), intent(in) :: prob !! Quantile probability in the closed interval [0, 1].
      real(dp) :: q
      real(dp), allocatable :: work(:)
      real(dp) :: h
      real(dp) :: frac
      integer :: j
      integer :: n

      n = size(x)
      if (n <= 0) then
         q = 0.0_dp
         return
      end if

      allocate(work(n))
      work = x
      call sort_real(work)
      if (n == 1) then
         q = work(1)
         return
      end if

      h = 1.0_dp + real(n - 1, dp) * max(0.0_dp, min(1.0_dp, prob))
      j = int(floor(h))
      frac = h - real(j, dp)
      if (j >= n) then
         q = work(n)
      else
         q = (1.0_dp - frac) * work(j) + frac * work(j + 1)
      end if
   end function quantile_type7

   pure subroutine bspline_design(knots, x, degree, deriv, values)
      real(dp), intent(in) :: knots(:) !! Nondecreasing complete B-spline knot sequence.
      real(dp), intent(in) :: x !! Evaluation point in the spline domain.
      integer, intent(in) :: degree !! Polynomial degree, at least zero.
      integer, intent(in) :: deriv !! Derivative order in [0, degree].
      real(dp), intent(out) :: values(:) !! B-spline basis or derivative values, length size(knots)-degree-1.
      real(dp), allocatable :: table(:, :, :)
      real(dp) :: xx
      real(dp) :: den_left
      real(dp) :: den_right
      integer :: d
      integer :: i
      integer :: m
      integer :: r
      integer :: ncoef
      integer :: nlev

      m = size(knots)
      ncoef = m - degree - 1
      values = 0.0_dp
      if (degree < 0 .or. deriv < 0 .or. deriv > degree) return
      if (ncoef <= 0 .or. size(values) /= ncoef) return

      xx = x
      if (xx >= knots(m)) xx = ieee_next_after(knots(m), -huge(xx))
      nlev = m - 1
      allocate(table(0:degree, 0:degree, nlev))
      table = 0.0_dp

      do i = 1, m - 1
         if (knots(i) <= xx .and. xx < knots(i + 1)) table(0, 0, i) = 1.0_dp
      end do

      do d = 1, degree
         do i = 1, m - d - 1
            den_left = knots(i + d) - knots(i)
            den_right = knots(i + d + 1) - knots(i + 1)
            if (den_left > 0.0_dp) then
               table(d, 0, i) = table(d, 0, i) + (xx - knots(i)) * table(d - 1, 0, i) / den_left
            end if
            if (den_right > 0.0_dp) then
               table(d, 0, i) = table(d, 0, i) + (knots(i + d + 1) - xx) * table(d - 1, 0, i + 1) / den_right
            end if
         end do
         do r = 1, d
            do i = 1, m - d - 1
               den_left = knots(i + d) - knots(i)
               den_right = knots(i + d + 1) - knots(i + 1)
               if (den_left > 0.0_dp) then
                  table(d, r, i) = table(d, r, i) + real(d, dp) * table(d - 1, r - 1, i) / den_left
               end if
               if (den_right > 0.0_dp) then
                  table(d, r, i) = table(d, r, i) - real(d, dp) * table(d - 1, r - 1, i + 1) / den_right
               end if
            end do
         end do
      end do

      values = table(degree, deriv, 1:ncoef)
   end subroutine bspline_design

end module pbs_splines
