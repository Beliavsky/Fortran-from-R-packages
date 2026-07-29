! SPDX-License-Identifier: GPL-3.0-only
!
! Copyright (C) 2017 Bernhard Pfaff
! Modern Fortran translation of the computational algorithms in mcrp.
module mcrp_moments
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use mcrp_kinds, only : dp, mcrp_success, mcrp_invalid_shape, &
      mcrp_invalid_argument
   implicit none
   private

   public :: m2, m3, m4
   public :: pm2, pm3, pm4
   public :: dm2, dm3, dm4
   public :: cm2, cm3, cm4
   public :: port_risk, port_risk_deriv, port_risk_contrib
   public :: port_skew, port_skew_deriv, port_skew_contrib
   public :: port_kurt, port_kurt_deriv, port_kurt_contrib

contains

   function m2(r, status) result(ans)
      real(dp), intent(in) :: r(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:, :)
      real(dp), allocatable :: rc(:, :)
      integer :: nobs, nasset

      nobs = size(r, 1)
      nasset = size(r, 2)
      if (nobs < 2 .or. nasset < 1) then
         allocate(ans(0, 0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if

      rc = centered_returns(r)
      ans = matmul(transpose(rc), rc) / real(nobs - 1, dp)
      call set_status(status, mcrp_success)
   end function m2

   function m3(r, status) result(ans)
      real(dp), intent(in) :: r(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:, :)
      real(dp), allocatable :: rc(:, :)
      integer :: i, a, b, c, col, nobs, nasset

      nobs = size(r, 1)
      nasset = size(r, 2)
      if (nobs < 1 .or. nasset < 1) then
         allocate(ans(0, 0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if

      rc = centered_returns(r)
      allocate(ans(nasset, nasset**2), source=0.0_dp)
      do i = 1, nobs
         do a = 1, nasset
            do b = 1, nasset
               do c = 1, nasset
                  col = (b - 1) * nasset + c
                  ans(a, col) = ans(a, col) + &
                     rc(i, a) * rc(i, b) * rc(i, c)
               end do
            end do
         end do
      end do
      ans = ans / real(nobs, dp)
      call set_status(status, mcrp_success)
   end function m3

   function m4(r, status) result(ans)
      real(dp), intent(in) :: r(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:, :)
      real(dp), allocatable :: rc(:, :)
      integer :: i, a, b, c, d, col, nobs, nasset

      nobs = size(r, 1)
      nasset = size(r, 2)
      if (nobs < 1 .or. nasset < 1) then
         allocate(ans(0, 0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if

      rc = centered_returns(r)
      allocate(ans(nasset, nasset**3), source=0.0_dp)
      do i = 1, nobs
         do a = 1, nasset
            do b = 1, nasset
               do c = 1, nasset
                  do d = 1, nasset
                     col = ((b - 1) * nasset + c - 1) * nasset + d
                     ans(a, col) = ans(a, col) + rc(i, a) * rc(i, b) * &
                        rc(i, c) * rc(i, d)
                  end do
               end do
            end do
         end do
      end do
      ans = ans / real(nobs, dp)
      call set_status(status, mcrp_success)
   end function m4

   function pm2(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans
      real(dp), allocatable :: me2(:, :)
      integer :: istat

      if (.not. weights_match(r, w)) then
         ans = quiet_nan()
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me2 = m2(r, istat)
      ans = dot_product(w, matmul(me2, w))
      call set_status(status, istat)
   end function pm2

   function pm3(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans
      real(dp), allocatable :: me3(:, :), ww(:)
      integer :: istat

      if (.not. weights_match(r, w)) then
         ans = quiet_nan()
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me3 = m3(r, istat)
      ww = kron2(w, w)
      ans = dot_product(w, matmul(me3, ww))
      call set_status(status, istat)
   end function pm3

   function pm4(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans
      real(dp), allocatable :: me4(:, :), www(:)
      integer :: istat

      if (.not. weights_match(r, w)) then
         ans = quiet_nan()
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me4 = m4(r, istat)
      www = kron3(w, w, w)
      ans = dot_product(w, matmul(me4, www))
      call set_status(status, istat)
   end function pm4

   function dm2(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: me2(:, :)
      integer :: istat

      if (.not. weights_match(r, w)) then
         allocate(ans(0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me2 = m2(r, istat)
      ans = 2.0_dp * matmul(me2, w)
      call set_status(status, istat)
   end function dm2

   function dm3(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: me3(:, :), ww(:)
      integer :: istat

      if (.not. weights_match(r, w)) then
         allocate(ans(0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me3 = m3(r, istat)
      ww = kron2(w, w)
      ans = 3.0_dp * matmul(me3, ww)
      call set_status(status, istat)
   end function dm3

   function dm4(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: me4(:, :), www(:)
      integer :: istat

      if (.not. weights_match(r, w)) then
         allocate(ans(0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      me4 = m4(r, istat)
      www = kron3(w, w, w)
      ans = 4.0_dp * matmul(me4, www)
      call set_status(status, istat)
   end function dm4

   function cm2(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      logical :: as_percentage
      real(dp) :: moment
      integer :: istat

      as_percentage = .true.
      if (present(percentage)) as_percentage = percentage
      ans = w * dm2(r, w, istat)
      if (istat /= mcrp_success) then
         call set_status(status, istat)
         return
      end if
      if (as_percentage) then
         moment = pm2(r, w)
         if (abs(moment) <= tiny(1.0_dp)) then
            ans = quiet_nan()
            call set_status(status, mcrp_invalid_argument)
            return
         end if
         ans = ans / moment
      end if
      ans = ans / 2.0_dp
      call set_status(status, mcrp_success)
   end function cm2

   function cm3(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      logical :: as_percentage
      real(dp) :: moment
      integer :: istat

      as_percentage = .true.
      if (present(percentage)) as_percentage = percentage
      ans = w * dm3(r, w, istat)
      if (istat /= mcrp_success) then
         call set_status(status, istat)
         return
      end if
      if (as_percentage) then
         moment = pm3(r, w)
         if (abs(moment) <= tiny(1.0_dp)) then
            ans = quiet_nan()
            call set_status(status, mcrp_invalid_argument)
            return
         end if
         ans = ans / moment
      end if
      ans = ans / 3.0_dp
      call set_status(status, mcrp_success)
   end function cm3

   function cm4(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      logical :: as_percentage
      real(dp) :: moment
      integer :: istat

      as_percentage = .true.
      if (present(percentage)) as_percentage = percentage
      ans = w * dm4(r, w, istat)
      if (istat /= mcrp_success) then
         call set_status(status, istat)
         return
      end if
      if (as_percentage) then
         moment = pm4(r, w)
         if (abs(moment) <= tiny(1.0_dp)) then
            ans = quiet_nan()
            call set_status(status, mcrp_invalid_argument)
            return
         end if
         ans = ans / moment
      end if
      ans = ans / 4.0_dp
      call set_status(status, mcrp_success)
   end function cm4

   function port_risk(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans
      ans = pm2(r, w, status)
   end function port_risk

   function port_risk_deriv(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      ans = dm2(r, w, status)
   end function port_risk_deriv

   function port_risk_contrib(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      ans = cm2(r, w, percentage, status)
   end function port_risk_contrib

   function port_skew(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans, var, third
      integer :: istat

      var = pm2(r, w, istat)
      if (istat /= mcrp_success .or. var <= 0.0_dp) then
         ans = quiet_nan()
         call set_status(status, mcrp_invalid_argument)
         return
      end if
      third = pm3(r, w)
      ans = third / var**1.5_dp
      call set_status(status, mcrp_success)
   end function port_skew

   function port_skew_deriv(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: d2(:), d3(:)
      real(dp) :: var, third
      integer :: istat

      if (.not. weights_match(r, w)) then
         allocate(ans(0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      var = pm2(r, w, istat)
      if (istat /= mcrp_success .or. var <= 0.0_dp) then
         allocate(ans(size(w)), source=quiet_nan())
         call set_status(status, mcrp_invalid_argument)
         return
      end if
      third = pm3(r, w)
      d2 = dm2(r, w)
      d3 = dm3(r, w)
      ans = (var**1.5_dp * d3 - third * sqrt(var) * d2) / var**3
      call set_status(status, mcrp_success)
   end function port_skew_deriv

   function port_skew_contrib(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      logical :: as_percentage
      real(dp) :: risk
      integer :: istat

      as_percentage = .true.
      if (present(percentage)) as_percentage = percentage
      ans = w * port_skew_deriv(r, w, istat)
      if (istat /= mcrp_success) then
         call set_status(status, istat)
         return
      end if
      if (as_percentage) then
         risk = port_skew(r, w)
         if (abs(risk) <= tiny(1.0_dp)) then
            ans = quiet_nan()
            call set_status(status, mcrp_invalid_argument)
            return
         end if
         ans = ans / risk
      end if
      call set_status(status, mcrp_success)
   end function port_skew_contrib

   function port_kurt(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp) :: ans, var, fourth
      integer :: istat

      var = pm2(r, w, istat)
      if (istat /= mcrp_success .or. var <= 0.0_dp) then
         ans = quiet_nan()
         call set_status(status, mcrp_invalid_argument)
         return
      end if
      fourth = pm4(r, w)
      ans = fourth / var**2
      call set_status(status, mcrp_success)
   end function port_kurt

   function port_kurt_deriv(r, w, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: d2(:), d4(:)
      real(dp) :: var, fourth
      integer :: istat

      if (.not. weights_match(r, w)) then
         allocate(ans(0))
         call set_status(status, mcrp_invalid_shape)
         return
      end if
      var = pm2(r, w, istat)
      if (istat /= mcrp_success .or. var <= 0.0_dp) then
         allocate(ans(size(w)), source=quiet_nan())
         call set_status(status, mcrp_invalid_argument)
         return
      end if
      fourth = pm4(r, w)
      d2 = dm2(r, w)
      d4 = dm4(r, w)
      ans = (var * d4 - fourth * d2) / (2.0_dp * var**3)
      call set_status(status, mcrp_success)
   end function port_kurt_deriv

   function port_kurt_contrib(r, w, percentage, status) result(ans)
      real(dp), intent(in) :: r(:, :), w(:)
      logical, intent(in), optional :: percentage
      integer, intent(out), optional :: status
      real(dp), allocatable :: ans(:)
      logical :: as_percentage
      real(dp) :: risk
      integer :: istat

      as_percentage = .true.
      if (present(percentage)) as_percentage = percentage
      ans = w * port_kurt_deriv(r, w, istat)
      if (istat /= mcrp_success) then
         call set_status(status, istat)
         return
      end if
      if (as_percentage) then
         risk = port_kurt(r, w)
         if (abs(risk) <= tiny(1.0_dp)) then
            ans = quiet_nan()
            call set_status(status, mcrp_invalid_argument)
            return
         end if
         ans = ans / risk
      end if
      call set_status(status, mcrp_success)
   end function port_kurt_contrib

   function centered_returns(r) result(rc)
      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: rc(:, :)
      real(dp), allocatable :: means(:)
      integer :: j

      means = sum(r, dim=1) / real(size(r, 1), dp)
      allocate(rc(size(r, 1), size(r, 2)))
      do j = 1, size(r, 2)
         rc(:, j) = r(:, j) - means(j)
      end do
   end function centered_returns

   function kron2(a, b) result(ans)
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable :: ans(:)
      integer :: i, j, k

      allocate(ans(size(a) * size(b)))
      k = 0
      do i = 1, size(a)
         do j = 1, size(b)
            k = k + 1
            ans(k) = a(i) * b(j)
         end do
      end do
   end function kron2

   function kron3(a, b, c) result(ans)
      real(dp), intent(in) :: a(:), b(:), c(:)
      real(dp), allocatable :: ans(:)
      integer :: i, j, k, l

      allocate(ans(size(a) * size(b) * size(c)))
      l = 0
      do i = 1, size(a)
         do j = 1, size(b)
            do k = 1, size(c)
               l = l + 1
               ans(l) = a(i) * b(j) * c(k)
            end do
         end do
      end do
   end function kron3

   logical function weights_match(r, w)
      real(dp), intent(in) :: r(:, :), w(:)
      weights_match = size(r, 1) >= 1 .and. size(r, 2) == size(w) .and. &
         size(w) >= 1
   end function weights_match

   real(dp) function quiet_nan()
      quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   subroutine set_status(status, value)
      integer, intent(out), optional :: status
      integer, intent(in) :: value
      if (present(status)) status = value
   end subroutine set_status

end module mcrp_moments
