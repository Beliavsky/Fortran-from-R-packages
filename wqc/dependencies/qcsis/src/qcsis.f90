module qcsis_mod
   use qcsis_kinds, only : dp
   use qcsis_statistics, only : descending_order, mean_value, quantile_type7_sorted, sample_sd, sort_real
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: qc_result
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: rho(:)
   end type qc_result

   type, public :: screening_result
      real(dp), allocatable :: w(:)
      integer, allocatable :: selected(:)
   end type screening_result

   interface qcsis
      module procedure qcsis_default
      module procedure qcsis_at_tau
   end interface qcsis

   public :: dp
   public :: qc
   public :: cqc
   public :: qcsis
   public :: cqcsis

contains

   ! Compute quantile correlations for one predictor at supplied probabilities.
   function qc(x, y, tau, stat, errmsg) result(fit)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: tau(:)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(qc_result) :: fit

      real(dp), allocatable :: sorted_y(:)
      real(dp) :: sx, xbar
      integer :: i, status
      character(len=:), allocatable :: message

      call initialize_status(stat, errmsg)
      call validate_vectors(x, y, status, message)
      if (status /= 0) then
         call fail_qc(fit, status, message, stat, errmsg)
         return
      end if

      call validate_tau(tau, status, message)
      if (status /= 0) then
         call fail_qc(fit, status, message, stat, errmsg)
         return
      end if

      sx = sample_sd(x)
      if (.not. ieee_is_finite(sx) .or. sx <= 0.0_dp) then
         call fail_qc(fit, 4, "x must have a finite, positive sample standard deviation", stat, errmsg)
         return
      end if

      allocate(fit%tau(size(tau)), fit%rho(size(tau)), sorted_y(size(y)))
      fit%tau = tau
      sorted_y = y
      call sort_real(sorted_y)
      xbar = mean_value(x)

      do i = 1, size(tau)
         fit%rho(i) = correlation_at_tau(x, y, sorted_y, xbar, sx, tau(i))
      end do
   end function qc


   ! Compute the mean quantile correlation over tau = 1/n, ..., (n-1)/n.
   function cqc(x, y, stat, errmsg) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      real(dp) :: value

      real(dp), allocatable :: sorted_y(:)
      real(dp) :: sx, tau, xbar
      integer :: i, n, status
      character(len=:), allocatable :: message

      call initialize_status(stat, errmsg)
      call validate_vectors(x, y, status, message)
      if (status /= 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         call set_status(status, message, stat, errmsg)
         return
      end if

      sx = sample_sd(x)
      if (.not. ieee_is_finite(sx) .or. sx <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         call set_status(4, "x must have a finite, positive sample standard deviation", stat, errmsg)
         return
      end if

      n = size(y)
      allocate(sorted_y(n))
      sorted_y = y
      call sort_real(sorted_y)
      xbar = mean_value(x)
      value = 0.0_dp
      do i = 1, n - 1
         tau = real(i, dp) / real(n, dp)
         value = value + correlation_at_tau(x, y, sorted_y, xbar, sx, tau)
      end do
      value = value / real(n - 1, dp)
   end function cqc


   ! Screen predictors using the default quantile grid from the R package.
   function qcsis_default(x, y, d, stat, errmsg) result(fit)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: d
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(screening_result) :: fit

      real(dp), allocatable :: tau(:)
      integer :: i, n

      n = size(y)
      allocate(tau(max(0, n - 1)))
      do i = 1, n - 1
         tau(i) = real(i, dp) / real(n, dp)
      end do
      fit = qcsis_at_tau(x, y, tau, d, stat, errmsg)
   end function qcsis_default


   ! Screen predictors using a caller-supplied quantile grid.
   function qcsis_at_tau(x, y, tau, d, stat, errmsg) result(fit)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: tau(:)
      integer, intent(in) :: d
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(screening_result) :: fit

      real(dp), allocatable :: psi(:), sorted_y(:), sx(:), xbar(:)
      real(dp) :: denominator, q, rho
      integer, allocatable :: order(:)
      integer :: i, j, n, p, status
      character(len=:), allocatable :: message

      call initialize_status(stat, errmsg)
      call validate_screening_inputs(x, y, d, status, message)
      if (status /= 0) then
         call fail_screening(fit, status, message, stat, errmsg)
         return
      end if

      call validate_tau(tau, status, message)
      if (status /= 0) then
         call fail_screening(fit, status, message, stat, errmsg)
         return
      end if

      n = size(x, 1)
      p = size(x, 2)
      allocate(fit%w(p), psi(n), sorted_y(n), sx(p), xbar(p))
      fit%w = 0.0_dp
      sorted_y = y
      call sort_real(sorted_y)

      do j = 1, p
         xbar(j) = mean_value(x(:, j))
         sx(j) = sample_sd(x(:, j))
         if (.not. ieee_is_finite(sx(j)) .or. sx(j) <= 0.0_dp) then
            call fail_screening(fit, 5, "every column of x must have a finite, positive sample standard deviation", &
               stat, errmsg)
            return
         end if
      end do

      do i = 1, size(tau)
         q = quantile_type7_sorted(sorted_y, tau(i))
         psi = tau(i)
         where (y < q) psi = psi - 1.0_dp
         denominator = real(n, dp) * sqrt(tau(i) - tau(i) * tau(i))
         do j = 1, p
            rho = centered_dot(x(:, j), xbar(j), psi) / (denominator * sx(j))
            fit%w(j) = fit%w(j) + rho * rho
         end do
      end do
      fit%w = fit%w / real(size(tau), dp)

      order = descending_order(fit%w)
      allocate(fit%selected(d))
      fit%selected = order(:d)
   end function qcsis_at_tau


   ! Screen predictors by absolute composite quantile correlation.
   function cqcsis(x, y, d, stat, errmsg) result(fit)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: d
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg
      type(screening_result) :: fit

      real(dp), allocatable :: psi(:), sorted_y(:), sx(:), xbar(:)
      real(dp) :: denominator, q, tau
      integer, allocatable :: order(:)
      integer :: i, j, n, p, status
      character(len=:), allocatable :: message

      call initialize_status(stat, errmsg)
      call validate_screening_inputs(x, y, d, status, message)
      if (status /= 0) then
         call fail_screening(fit, status, message, stat, errmsg)
         return
      end if

      n = size(x, 1)
      p = size(x, 2)
      allocate(fit%w(p), psi(n), sorted_y(n), sx(p), xbar(p))
      fit%w = 0.0_dp
      sorted_y = y
      call sort_real(sorted_y)

      do j = 1, p
         xbar(j) = mean_value(x(:, j))
         sx(j) = sample_sd(x(:, j))
         if (.not. ieee_is_finite(sx(j)) .or. sx(j) <= 0.0_dp) then
            call fail_screening(fit, 5, "every column of x must have a finite, positive sample standard deviation", &
               stat, errmsg)
            return
         end if
      end do

      do i = 1, n - 1
         tau = real(i, dp) / real(n, dp)
         q = quantile_type7_sorted(sorted_y, tau)
         psi = tau
         where (y < q) psi = psi - 1.0_dp
         denominator = real(n, dp) * sqrt(tau - tau * tau)
         do j = 1, p
            fit%w(j) = fit%w(j) + centered_dot(x(:, j), xbar(j), psi) / (denominator * sx(j))
         end do
      end do
      fit%w = abs(fit%w / real(n - 1, dp))

      order = descending_order(fit%w)
      allocate(fit%selected(d))
      fit%selected = order(:d)
   end function cqcsis


   function correlation_at_tau(x, y, sorted_y, xbar, sx, tau) result(rho)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in) :: sorted_y(:)
      real(dp), intent(in) :: xbar
      real(dp), intent(in) :: sx
      real(dp), intent(in) :: tau
      real(dp) :: rho

      real(dp) :: q, psi
      integer :: i, n

      n = size(y)
      q = quantile_type7_sorted(sorted_y, tau)
      rho = 0.0_dp
      do i = 1, n
         psi = tau
         if (y(i) < q) psi = psi - 1.0_dp
         rho = rho + (x(i) - xbar) * psi
      end do
      rho = rho / (real(n, dp) * sqrt(tau - tau * tau) * sx)
   end function correlation_at_tau


   pure function centered_dot(x, xbar, weights) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: xbar
      real(dp), intent(in) :: weights(:)
      real(dp) :: value

      integer :: i

      value = 0.0_dp
      do i = 1, size(x)
         value = value + (x(i) - xbar) * weights(i)
      end do
   end function centered_dot


   subroutine validate_vectors(x, y, stat, errmsg)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      stat = 0
      errmsg = ""
      if (size(x) /= size(y)) then
         stat = 1
         errmsg = "x and y must have the same length"
      else if (size(y) < 2) then
         stat = 2
         errmsg = "x and y must contain at least two observations"
      else if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         stat = 3
         errmsg = "x and y must contain only finite values"
      end if
   end subroutine validate_vectors


   subroutine validate_tau(tau, stat, errmsg)
      real(dp), intent(in) :: tau(:)
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      stat = 0
      errmsg = ""
      if (size(tau) == 0) then
         stat = 6
         errmsg = "tau must contain at least one probability"
      else if (.not. all(ieee_is_finite(tau))) then
         stat = 7
         errmsg = "tau must contain only finite probabilities"
      else if (any(tau <= 0.0_dp) .or. any(tau >= 1.0_dp)) then
         stat = 8
         errmsg = "every tau value must be strictly between zero and one"
      end if
   end subroutine validate_tau


   subroutine validate_screening_inputs(x, y, d, stat, errmsg)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: d
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      stat = 0
      errmsg = ""
      if (size(x, 1) /= size(y)) then
         stat = 1
         errmsg = "the number of rows in x must equal the length of y"
      else if (size(y) < 2) then
         stat = 2
         errmsg = "x and y must contain at least two observations"
      else if (size(x, 2) < 1) then
         stat = 9
         errmsg = "x must contain at least one predictor column"
      else if (d < 1 .or. d > size(x, 2)) then
         stat = 10
         errmsg = "d must be between one and the number of predictor columns"
      else if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         stat = 3
         errmsg = "x and y must contain only finite values"
      end if
   end subroutine validate_screening_inputs


   subroutine initialize_status(stat, errmsg)
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      if (present(stat)) stat = 0
      if (present(errmsg)) errmsg = ""
   end subroutine initialize_status


   subroutine set_status(code, message, stat, errmsg)
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      if (present(stat)) stat = code
      if (present(errmsg)) errmsg = message
   end subroutine set_status


   subroutine fail_qc(fit, code, message, stat, errmsg)
      type(qc_result), intent(out) :: fit
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      allocate(fit%tau(0), fit%rho(0))
      call set_status(code, message, stat, errmsg)
   end subroutine fail_qc


   subroutine fail_screening(fit, code, message, stat, errmsg)
      type(screening_result), intent(out) :: fit
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      integer, optional, intent(out) :: stat
      character(len=:), allocatable, optional, intent(out) :: errmsg

      if (allocated(fit%w)) deallocate(fit%w)
      if (allocated(fit%selected)) deallocate(fit%selected)
      allocate(fit%w(0), fit%selected(0))
      call set_status(code, message, stat, errmsg)
   end subroutine fail_screening

end module qcsis_mod
