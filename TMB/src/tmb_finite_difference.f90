module tmb_finite_difference
   use tmb_kinds, only: dp
   implicit none
   private
   public :: gradient_fd, hessian_fd
   abstract interface
      pure function scalar_objective(x) result(value)
         import dp
         real(dp), intent(in) :: x(:) !! Parameter vector at which the objective is evaluated.
         real(dp) :: value
      end function scalar_objective
   end interface
contains
   pure subroutine gradient_fd(f, x, gradient, step)
      procedure(scalar_objective) :: f !! Pure scalar objective function.
      real(dp), intent(in) :: x(:) !! Parameter vector at which to approximate the gradient.
      real(dp), intent(out) :: gradient(:) !! Central-difference gradient, same length as x.
      real(dp), intent(in), optional :: step !! Positive relative finite-difference step; default sqrt(epsilon).
      real(dp) :: xp(size(x)), xm(size(x)), h, base_step
      integer :: i
      base_step = sqrt(epsilon(1.0_dp))
      if (present(step)) base_step = step
      do i = 1, size(x)
         h = base_step * max(1.0_dp, abs(x(i)))
         xp = x
         xm = x
         xp(i) = xp(i) + h
         xm(i) = xm(i) - h
         gradient(i) = (f(xp) - f(xm)) / (2.0_dp * h)
      end do
   end subroutine gradient_fd

   pure subroutine hessian_fd(f, x, hessian, step)
      procedure(scalar_objective) :: f !! Pure scalar objective function.
      real(dp), intent(in) :: x(:) !! Parameter vector at which to approximate the Hessian.
      real(dp), intent(out) :: hessian(:, :) !! Symmetric Hessian matrix with shape (size(x),size(x)).
      real(dp), intent(in), optional :: step !! Positive relative finite-difference step; default epsilon**0.25.
      real(dp) :: xpp(size(x)), xpm(size(x)), xmp(size(x)), xmm(size(x))
      real(dp) :: hi, hj, base_step, f0
      integer :: i, j
      base_step = epsilon(1.0_dp)**0.25_dp
      if (present(step)) base_step = step
      f0 = f(x)
      hessian = 0.0_dp
      do i = 1, size(x)
         hi = base_step * max(1.0_dp, abs(x(i)))
         xpp = x
         xpm = x
         xpp(i) = xpp(i) + hi
         xpm(i) = xpm(i) - hi
         hessian(i, i) = (f(xpp) - 2.0_dp * f0 + f(xpm)) / (hi * hi)
         do j = i + 1, size(x)
            hj = base_step * max(1.0_dp, abs(x(j)))
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i) + hi
            xpp(j) = xpp(j) + hj
            xpm(i) = xpm(i) + hi
            xpm(j) = xpm(j) - hj
            xmp(i) = xmp(i) - hi
            xmp(j) = xmp(j) + hj
            xmm(i) = xmm(i) - hi
            xmm(j) = xmm(j) - hj
            hessian(i, j) = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4.0_dp * hi * hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine hessian_fd
end module tmb_finite_difference
