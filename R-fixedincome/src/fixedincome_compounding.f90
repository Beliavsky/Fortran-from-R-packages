! SPDX-License-Identifier: MIT
module fixedincome_compounding
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use fixedincome_kinds, only : dp
   use fixedincome_types, only : COMPOUND_SIMPLE, COMPOUND_DISCRETE, COMPOUND_CONTINUOUS, &
      FI_OK, FI_INVALID_ARGUMENT
   implicit none
   private

   public :: compounding, compounding_name
   public :: compound, implied_rate, discount

   interface compound
      module procedure compound_scalar
      module procedure compound_vector
   end interface compound

   interface implied_rate
      module procedure implied_rate_scalar
      module procedure implied_rate_vector
   end interface implied_rate

   interface discount
      module procedure discount_scalar
      module procedure discount_vector
   end interface discount

contains

   pure integer function compounding(name) result(method)
      character(len=*), intent(in) :: name
      select case (lower(trim(adjustl(name))))
      case ('simple')
         method = COMPOUND_SIMPLE
      case ('discrete')
         method = COMPOUND_DISCRETE
      case ('continuous')
         method = COMPOUND_CONTINUOUS
      case default
         method = 0
      end select
   end function compounding

   pure function compounding_name(method) result(name)
      integer, intent(in) :: method
      character(len=16) :: name
      select case (method)
      case (COMPOUND_SIMPLE)
         name = 'simple'
      case (COMPOUND_DISCRETE)
         name = 'discrete'
      case (COMPOUND_CONTINUOUS)
         name = 'continuous'
      case default
         name = 'unknown'
      end select
   end function compounding_name

   function compound_scalar(method, time, rate, status) result(factor)
      integer, intent(in) :: method
      real(dp), intent(in) :: time, rate
      integer, intent(out), optional :: status
      real(dp) :: factor
      factor = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(rate)) then
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      select case (method)
      case (COMPOUND_SIMPLE)
         factor = 1.0_dp + rate * time
      case (COMPOUND_DISCRETE)
         if (1.0_dp + rate < 0.0_dp .and. abs(time - anint(time)) > 32.0_dp * epsilon(time)) then
            if (present(status)) status = FI_INVALID_ARGUMENT
            return
         end if
         factor = (1.0_dp + rate)**time
      case (COMPOUND_CONTINUOUS)
         factor = exp(rate * time)
      case default
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = FI_OK
   end function compound_scalar

   function compound_vector(method, time, rate, status) result(factor)
      integer, intent(in) :: method
      real(dp), intent(in) :: time(:), rate(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: factor(:)
      integer :: i, stat_i
      if (size(time) /= size(rate)) then
         allocate(factor(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      allocate(factor(size(time)))
      do i = 1, size(time)
         factor(i) = compound_scalar(method, time(i), rate(i), stat_i)
         if (stat_i /= FI_OK) then
            if (present(status)) status = stat_i
            return
         end if
      end do
      if (present(status)) status = FI_OK
   end function compound_vector

   function implied_rate_scalar(method, time, factor, status) result(rate)
      integer, intent(in) :: method
      real(dp), intent(in) :: time, factor
      integer, intent(out), optional :: status
      real(dp) :: rate
      rate = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(factor)) then
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      if (abs(time) <= tiny(1.0_dp)) then
         if (abs(factor - 1.0_dp) <= 32.0_dp * epsilon(factor)) then
            rate = 0.0_dp
            if (present(status)) status = FI_OK
         else
            if (present(status)) status = FI_INVALID_ARGUMENT
         end if
         return
      end if
      select case (method)
      case (COMPOUND_SIMPLE)
         rate = (factor - 1.0_dp) / time
      case (COMPOUND_DISCRETE)
         if (factor < 0.0_dp) then
            if (present(status)) status = FI_INVALID_ARGUMENT
            return
         end if
         rate = factor**(1.0_dp / time) - 1.0_dp
      case (COMPOUND_CONTINUOUS)
         if (factor <= 0.0_dp) then
            if (present(status)) status = FI_INVALID_ARGUMENT
            return
         end if
         rate = log(factor) / time
      case default
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = FI_OK
   end function implied_rate_scalar

   function implied_rate_vector(method, time, factor, status) result(rate)
      integer, intent(in) :: method
      real(dp), intent(in) :: time(:), factor(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: rate(:)
      integer :: i, stat_i
      if (size(time) /= size(factor)) then
         allocate(rate(0))
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if
      allocate(rate(size(time)))
      do i = 1, size(time)
         rate(i) = implied_rate_scalar(method, time(i), factor(i), stat_i)
         if (stat_i /= FI_OK) then
            if (present(status)) status = stat_i
            return
         end if
      end do
      if (present(status)) status = FI_OK
   end function implied_rate_vector

   function discount_scalar(method, time, rate, status) result(value)
      integer, intent(in) :: method
      real(dp), intent(in) :: time, rate
      integer, intent(out), optional :: status
      real(dp) :: value, factor
      integer :: stat_i
      factor = compound_scalar(method, time, rate, stat_i)
      if (stat_i /= FI_OK .or. abs(factor) <= tiny(1.0_dp)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
      else
         value = 1.0_dp / factor
         if (present(status)) status = FI_OK
      end if
   end function discount_scalar

   function discount_vector(method, time, rate, status) result(value)
      integer, intent(in) :: method
      real(dp), intent(in) :: time(:), rate(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: factor(:)
      integer :: stat_i
      factor = compound_vector(method, time, rate, stat_i)
      if (stat_i /= FI_OK .or. any(abs(factor) <= tiny(1.0_dp))) then
         allocate(value(size(factor)))
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
      else
         allocate(value(size(factor)))
         value = 1.0_dp / factor
         if (present(status)) status = FI_OK
      end if
   end function discount_vector

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, c
      do i = 1, len(text)
         c = iachar(text(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            out(i:i) = achar(c + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower

end module fixedincome_compounding
