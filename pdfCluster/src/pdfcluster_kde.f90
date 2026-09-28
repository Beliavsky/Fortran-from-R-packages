module pdfcluster_kde
   use pdfcluster_kinds, only : dp
   use pdfcluster_types, only : kepdf_result
   use pdfcluster_stats, only : sample_sd_columns
   implicit none
   private

   public :: h_norm, hprop2f, kepdf_estimate, kepdf_fit

contains

   pure subroutine h_norm(x, h)
      real(dp), intent(in) :: x(:, :) !! Numeric sample, with observations in rows and dimensions in columns.
      real(dp), allocatable, intent(out) :: h(:) !! Normal-reference bandwidth for each data dimension.
      real(dp), allocatable :: sd(:)
      real(dp) :: factor
      integer :: d
      integer :: n

      n = size(x, 1)
      d = size(x, 2)
      allocate(h(d), sd(d))
      if (n <= 0 .or. d <= 0) then
         h = 0.0_dp
         return
      end if
      call sample_sd_columns(x, sd)
      factor = (4.0_dp/(real(d + 2, dp)*real(n, dp)))**(1.0_dp/real(d + 4, dp))
      h = sd*factor
   end subroutine h_norm

   pure subroutine kepdf_estimate(x, eval_points, h, hx, kernel, estimate, success)
      real(dp), intent(in) :: x(:, :) !! Kernel centers, one observation per row.
      real(dp), intent(in) :: eval_points(:, :) !! Evaluation locations with the same number of columns as x.
      real(dp), intent(in), optional :: h(:) !! Fixed diagonal bandwidth; omitted when row-specific hx is supplied.
      real(dp), intent(in), optional :: hx(:, :) !! Row-specific diagonal bandwidth matrix, shaped like x.
      character(len=*), intent(in), optional :: kernel !! Kernel name: "gaussian" or "t7"; default is Gaussian.
      real(dp), allocatable, intent(out) :: estimate(:) !! Estimated density at each evaluation point.
      logical, intent(out) :: success !! True when dimensions and strictly positive bandwidths are valid.
      character(len=12) :: kern
      real(dp), allocatable :: hm(:, :)
      real(dp) :: cons
      real(dp) :: prod_inv_h
      real(dp) :: prod_term
      real(dp) :: s
      real(dp) :: u
      integer :: d
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: ne

      success = .false.
      n = size(x, 1)
      d = size(x, 2)
      ne = size(eval_points, 1)
      allocate(estimate(ne))
      estimate = 0.0_dp
      if (n <= 0 .or. d <= 0 .or. size(eval_points, 2) /= d) return
      allocate(hm(n, d))
      if (present(hx)) then
         if (size(hx, 1) /= n .or. size(hx, 2) /= d) return
         hm = hx
      else if (present(h)) then
         if (size(h) /= d) return
         do k = 1, n
            hm(k, :) = h
         end do
      else
         return
      end if
      if (any(hm <= 0.0_dp)) return
      kern = "gaussian"
      if (present(kernel)) kern = adjustl(kernel)

      select case (trim(kern))
      case ("gaussian")
         cons = 0.3989423_dp**d/real(n, dp)
         do i = 1, ne
            s = 0.0_dp
            do k = 1, n
               prod_inv_h = product(1.0_dp/hm(k, :))
               prod_term = 0.0_dp
               do j = 1, d
                  u = (x(k, j) - eval_points(i, j))/(hm(k, j)*sqrt(2.0_dp))
                  prod_term = prod_term + u*u
               end do
               s = s + exp(-prod_term)*prod_inv_h
            end do
            estimate(i) = s*cons
         end do
      case ("t7")
         cons = 0.3849915_dp**d/real(n, dp)
         do i = 1, ne
            s = 0.0_dp
            do k = 1, n
               prod_inv_h = product(1.0_dp/hm(k, :))
               prod_term = 1.0_dp
               do j = 1, d
                  u = (x(k, j) - eval_points(i, j))/hm(k, j)
                  prod_term = prod_term*(7.0_dp + u*u)
               end do
               s = s + 7.0_dp**(4*d)*prod_inv_h/(prod_term**4)
            end do
            estimate(i) = s*cons
         end do
      case default
         return
      end select
      success = .true.
   end subroutine kepdf_estimate

   pure subroutine hprop2f(x, h, hx, alpha, kernel, success)
      real(dp), intent(in) :: x(:, :) !! Sample used for the pilot density and adaptive bandwidth calculation.
      real(dp), intent(in) :: h(:) !! Positive fixed pilot bandwidth vector.
      real(dp), allocatable, intent(out) :: hx(:, :) !! Adaptive bandwidth matrix, one row per observation.
      real(dp), intent(in), optional :: alpha !! Abramson exponent; default is 0.5.
      character(len=*), intent(in), optional :: kernel !! Kernel name used for the pilot estimate.
      logical, intent(out) :: success !! True when the pilot density and adaptive bandwidths are valid.
      real(dp), allocatable :: pilot(:)
      real(dp) :: a
      real(dp) :: g
      real(dp) :: lambda
      integer :: i
      integer :: n

      success = .false.
      n = size(x, 1)
      allocate(hx(n, size(x, 2)))
      if (n <= 0 .or. size(h) /= size(x, 2) .or. any(h <= 0.0_dp)) return
      a = 0.5_dp
      if (present(alpha)) a = alpha
      call kepdf_estimate(x, x, h=h, kernel=kernel, estimate=pilot, success=success)
      if (.not. success) return
      if (any(pilot <= 0.0_dp)) then
         success = .false.
         return
      end if
      g = exp(sum(log(pilot))/real(n, dp))
      do i = 1, n
         lambda = (pilot(i)/g)**(-a)
         hx(i, :) = lambda*h
      end do
      success = .true.
   end subroutine hprop2f

   pure subroutine kepdf_fit(x, eval_points, result, h, hx, kernel, bwtype, alpha, success)
      real(dp), intent(in) :: x(:, :) !! Sample used as kernel centers.
      real(dp), intent(in) :: eval_points(:, :) !! Locations where the density is evaluated.
      type(kepdf_result), intent(out) :: result !! Density estimate and bandwidth metadata.
      real(dp), intent(in), optional :: h(:) !! Fixed or pilot bandwidth vector; normal-reference bandwidth is the default.
      real(dp), intent(in), optional :: hx(:, :) !! Row-specific adaptive bandwidth matrix when already available.
      character(len=*), intent(in), optional :: kernel !! Kernel name, "gaussian" or "t7".
      character(len=*), intent(in), optional :: bwtype !! Bandwidth mode, "fixed" or "adaptive".
      real(dp), intent(in), optional :: alpha !! Adaptive-bandwidth exponent; default is 0.5.
      logical, intent(out) :: success !! True when the requested density estimate is computed.
      real(dp), allocatable :: h0(:)
      character(len=12) :: bw
      character(len=12) :: kern
      real(dp) :: a

      success = .false.
      kern = "gaussian"
      bw = "fixed"
      a = 0.5_dp
      if (present(kernel)) kern = adjustl(kernel)
      if (present(bwtype)) bw = adjustl(bwtype)
      if (present(alpha)) a = alpha
      if (present(h)) then
         allocate(h0(size(h)))
         h0 = h
      else
         call h_norm(x, h0)
      end if
      result%kernel = kern
      result%bwtype = bw
      result%alpha = a
      result%h = h0
      if (trim(bw) == "adaptive") then
         if (present(hx)) then
            allocate(result%hx(size(hx, 1), size(hx, 2)))
            result%hx = hx
         else
            call hprop2f(x, h0, result%hx, alpha=a, kernel=kern, success=success)
            if (.not. success) return
         end if
         call kepdf_estimate(x, eval_points, hx=result%hx, kernel=kern, estimate=result%estimate, success=success)
      else if (trim(bw) == "fixed") then
         allocate(result%hx(0, 0))
         call kepdf_estimate(x, eval_points, h=h0, kernel=kern, estimate=result%estimate, success=success)
      end if
   end subroutine kepdf_fit

end module pdfcluster_kde
