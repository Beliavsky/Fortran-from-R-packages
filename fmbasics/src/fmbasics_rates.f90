! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_rates
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT, FM_SIZE_MISMATCH, FM_DOMAIN_ERROR
   use fmbasics_dates, only : year_frac, make_date
   implicit none
   private

   real(dp), parameter, public :: COMPOUND_CONTINUOUS = huge(1.0_dp)

   type, public :: interest_rate_t
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: compounding(:)
      character(len=12), allocatable :: day_basis(:)
   contains
      procedure :: size => rate_size
   end type interest_rate_t

   type, public :: discount_factor_t
      real(dp), allocatable :: value(:)
      integer, allocatable :: start_date(:)
      integer, allocatable :: end_date(:)
   contains
      procedure :: size => df_size
   end type discount_factor_t

   public :: interest_rate, discount_factor
   public :: as_discount_factor, as_interest_rate
   public :: convert_interest_rate, compound_factor, implied_rate
   public :: is_valid_compounding
   public :: rate_add, rate_subtract, rate_multiply, rate_divide
   public :: discount_multiply, discount_divide

   interface interest_rate
      module procedure interest_rate_scalar
      module procedure interest_rate_homogeneous
      module procedure interest_rate_full
   end interface interest_rate

   interface discount_factor
      module procedure discount_factor_scalar
      module procedure discount_factor_vector
   end interface discount_factor

   interface as_discount_factor
      module procedure as_discount_factor_scalar_dates
      module procedure as_discount_factor_vector_dates
   end interface as_discount_factor

   interface as_interest_rate
      module procedure as_interest_rate_scalar_basis
      module procedure as_interest_rate_vector_basis
   end interface as_interest_rate

contains

   function interest_rate_scalar(value, compounding, day_basis, status) result(rate)
      real(dp), intent(in) :: value, compounding
      character(len=*), intent(in) :: day_basis
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      allocate(rate%value(1), rate%compounding(1), rate%day_basis(1))
      rate%value = value
      rate%compounding = compounding
      rate%day_basis = lower12(day_basis)
      call validate_rate(rate, status)
   end function interest_rate_scalar

   function interest_rate_homogeneous(value, compounding, day_basis, status) result(rate)
      real(dp), intent(in) :: value(:), compounding
      character(len=*), intent(in) :: day_basis
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      allocate(rate%value(size(value)), rate%compounding(size(value)), rate%day_basis(size(value)))
      rate%value = value
      rate%compounding = compounding
      rate%day_basis = lower12(day_basis)
      call validate_rate(rate, status)
   end function interest_rate_homogeneous

   function interest_rate_full(value, compounding, day_basis, status) result(rate)
      real(dp), intent(in) :: value(:), compounding(:)
      character(len=*), intent(in) :: day_basis(:)
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      integer :: i, n
      n = max(size(value), size(compounding), size(day_basis))
      allocate(rate%value(n), rate%compounding(n), rate%day_basis(n))
      do i = 1, n
         rate%value(i) = value(mod(i-1,size(value))+1)
         rate%compounding(i) = compounding(mod(i-1,size(compounding))+1)
         rate%day_basis(i) = lower12(day_basis(mod(i-1,size(day_basis))+1))
      end do
      call validate_rate(rate, status)
   end function interest_rate_full

   function discount_factor_scalar(value, d1, d2, status) result(df)
      real(dp), intent(in) :: value
      integer, intent(in) :: d1, d2
      integer, intent(out), optional :: status
      type(discount_factor_t) :: df
      allocate(df%value(1), df%start_date(1), df%end_date(1))
      df%value = value
      df%start_date = d1
      df%end_date = d2
      call validate_df(df, status)
   end function discount_factor_scalar

   function discount_factor_vector(value, d1, d2, status) result(df)
      real(dp), intent(in) :: value(:)
      integer, intent(in) :: d1(:), d2(:)
      integer, intent(out), optional :: status
      type(discount_factor_t) :: df
      integer :: i, n
      n = max(size(value), size(d1), size(d2))
      allocate(df%value(n), df%start_date(n), df%end_date(n))
      do i = 1, n
         df%value(i) = value(mod(i-1,size(value))+1)
         df%start_date(i) = d1(mod(i-1,size(d1))+1)
         df%end_date(i) = d2(mod(i-1,size(d2))+1)
      end do
      call validate_df(df, status)
   end function discount_factor_vector

   function as_discount_factor_scalar_dates(rate, d1, d2, status) result(df)
      type(interest_rate_t), intent(in) :: rate
      integer, intent(in) :: d1, d2
      integer, intent(out), optional :: status
      type(discount_factor_t) :: df
      integer, allocatable :: starts(:), ends(:)
      allocate(starts(rate%size()), ends(rate%size()))
      starts = d1
      ends = d2
      df = as_discount_factor_vector_dates(rate, starts, ends, status)
   end function as_discount_factor_scalar_dates

   function as_discount_factor_vector_dates(rate, d1, d2, status) result(df)
      type(interest_rate_t), intent(in) :: rate
      integer, intent(in) :: d1(:), d2(:)
      integer, intent(out), optional :: status
      type(discount_factor_t) :: df
      integer :: i, n, stat_i
      real(dp) :: term
      n = max(rate%size(), size(d1), size(d2))
      allocate(df%value(n), df%start_date(n), df%end_date(n))
      stat_i = FM_OK
      do i = 1, n
         df%start_date(i) = d1(mod(i-1,size(d1))+1)
         df%end_date(i) = d2(mod(i-1,size(d2))+1)
         term = year_frac(df%start_date(i), df%end_date(i), &
            rate%day_basis(mod(i-1,rate%size())+1), stat_i)
         if (stat_i /= FM_OK) exit
         df%value(i) = 1.0_dp / compound_factor( &
            rate%value(mod(i-1,rate%size())+1), &
            rate%compounding(mod(i-1,rate%size())+1), term, stat_i)
         if (stat_i /= FM_OK) exit
      end do
      if (stat_i == FM_OK) call validate_df(df, stat_i)
      if (present(status)) status = stat_i
   end function as_discount_factor_vector_dates

   function as_interest_rate_scalar_basis(df, compounding, day_basis, status) result(rate)
      type(discount_factor_t), intent(in) :: df
      real(dp), intent(in) :: compounding
      character(len=*), intent(in) :: day_basis
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      real(dp), allocatable :: compounds(:)
      character(len=12), allocatable :: bases(:)
      allocate(compounds(df%size()), bases(df%size()))
      compounds = compounding
      bases = lower12(day_basis)
      rate = as_interest_rate_vector_basis(df, compounds, bases, status)
   end function as_interest_rate_scalar_basis

   function as_interest_rate_vector_basis(df, compounding, day_basis, status) result(rate)
      type(discount_factor_t), intent(in) :: df
      real(dp), intent(in) :: compounding(:)
      character(len=*), intent(in) :: day_basis(:)
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      integer :: i, n, stat_i
      real(dp) :: term
      n = max(df%size(), size(compounding), size(day_basis))
      allocate(rate%value(n), rate%compounding(n), rate%day_basis(n))
      stat_i = FM_OK
      do i = 1, n
         rate%compounding(i) = compounding(mod(i-1,size(compounding))+1)
         rate%day_basis(i) = lower12(day_basis(mod(i-1,size(day_basis))+1))
         term = year_frac(df%start_date(mod(i-1,df%size())+1), &
            df%end_date(mod(i-1,df%size())+1), rate%day_basis(i), stat_i)
         rate%value(i) = implied_rate(df%value(mod(i-1,df%size())+1), &
            rate%compounding(i), term, stat_i)
         if (stat_i /= FM_OK) exit
      end do
      if (stat_i == FM_OK) call validate_rate(rate, stat_i)
      if (present(status)) status = stat_i
   end function as_interest_rate_vector_basis

   function convert_interest_rate(rate, compounding, day_basis, status) result(converted)
      type(interest_rate_t), intent(in) :: rate
      real(dp), intent(in), optional :: compounding
      character(len=*), intent(in), optional :: day_basis
      integer, intent(out), optional :: status
      type(interest_rate_t) :: converted
      type(discount_factor_t) :: df
      real(dp), allocatable :: compounds(:)
      character(len=12), allocatable :: bases(:)
      integer :: d1, d2
      d1 = make_date(2013, 1, 1)
      d2 = make_date(2014, 1, 1)
      allocate(compounds(rate%size()), bases(rate%size()))
      compounds = rate%compounding
      bases = rate%day_basis
      if (present(compounding)) compounds = compounding
      if (present(day_basis)) bases = lower12(day_basis)
      df = as_discount_factor(rate, d1, d2)
      converted = as_interest_rate(df, compounds, bases, status)
   end function convert_interest_rate

   real(dp) function compound_factor(rate, compounding, term, status) result(value)
      real(dp), intent(in) :: rate, compounding, term
      integer, intent(out), optional :: status
      if (.not. is_valid_compounding(compounding)) then
         value = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      if (is_continuous(compounding)) then
         value = exp(rate * term)
      else if (abs(compounding) < 1.0e-12_dp) then
         value = 1.0_dp + rate * term
      else if (abs(compounding + 1.0_dp) < 1.0e-12_dp) then
         value = 1.0_dp / (1.0_dp - rate * term)
      else
         if (1.0_dp + rate / compounding <= 0.0_dp) then
            value = 0.0_dp
            if (present(status)) status = FM_DOMAIN_ERROR
            return
         end if
         value = (1.0_dp + rate / compounding) ** (compounding * term)
      end if
      if (value <= 0.0_dp) then
         if (present(status)) status = FM_DOMAIN_ERROR
      else
         if (present(status)) status = FM_OK
      end if
   end function compound_factor

   real(dp) function implied_rate(discount, compounding, term, status) result(value)
      real(dp), intent(in) :: discount, compounding, term
      integer, intent(out), optional :: status
      if (discount <= 0.0_dp .or. abs(term) <= tiny(1.0_dp) .or. &
          .not. is_valid_compounding(compounding)) then
         value = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      if (is_continuous(compounding)) then
         value = -log(discount) / term
      else if (abs(compounding) < 1.0e-12_dp) then
         value = (1.0_dp / discount - 1.0_dp) / term
      else if (abs(compounding + 1.0_dp) < 1.0e-12_dp) then
         value = (1.0_dp - discount) / term
      else
         value = compounding * ((1.0_dp / discount) ** &
            (1.0_dp / (compounding * term)) - 1.0_dp)
      end if
      if (present(status)) status = FM_OK
   end function implied_rate

   pure logical function is_valid_compounding(compounding) result(value)
      real(dp), intent(in) :: compounding
      real(dp), parameter :: allowed(11) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
         4.0_dp, 6.0_dp, 12.0_dp, 24.0_dp, 52.0_dp, 365.0_dp]
      value = is_continuous(compounding) .or. any(abs(compounding - allowed) < 1.0e-12_dp)
   end function is_valid_compounding

   function rate_add(a, b, status) result(c)
      type(interest_rate_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(interest_rate_t) :: c
      c = rate_binary(a, b, 1, status)
   end function rate_add

   function rate_subtract(a, b, status) result(c)
      type(interest_rate_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(interest_rate_t) :: c
      c = rate_binary(a, b, 2, status)
   end function rate_subtract

   function rate_multiply(a, b, status) result(c)
      type(interest_rate_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(interest_rate_t) :: c
      c = rate_binary(a, b, 3, status)
   end function rate_multiply

   function rate_divide(a, b, status) result(c)
      type(interest_rate_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(interest_rate_t) :: c
      c = rate_binary(a, b, 4, status)
   end function rate_divide

   function rate_binary(a, b, operation, status) result(c)
      type(interest_rate_t), intent(in) :: a, b
      integer, intent(in) :: operation
      integer, intent(out), optional :: status
      type(interest_rate_t) :: c
      type(interest_rate_t) :: bi, bj
      integer :: i, n, stat_i
      n = max(a%size(), b%size())
      allocate(c%value(n), c%compounding(n), c%day_basis(n))
      stat_i = FM_OK
      do i = 1, n
         bi = interest_rate(b%value(mod(i-1,b%size())+1), &
            b%compounding(mod(i-1,b%size())+1), b%day_basis(mod(i-1,b%size())+1))
         bj = convert_interest_rate(bi, a%compounding(mod(i-1,a%size())+1), &
            a%day_basis(mod(i-1,a%size())+1), stat_i)
         c%compounding(i) = a%compounding(mod(i-1,a%size())+1)
         c%day_basis(i) = a%day_basis(mod(i-1,a%size())+1)
         select case (operation)
         case (1)
            c%value(i) = a%value(mod(i-1,a%size())+1) + bj%value(1)
         case (2)
            c%value(i) = a%value(mod(i-1,a%size())+1) - bj%value(1)
         case (3)
            c%value(i) = a%value(mod(i-1,a%size())+1) * bj%value(1)
         case (4)
            if (abs(bj%value(1)) <= tiny(1.0_dp)) then
               stat_i = FM_DOMAIN_ERROR
               c%value(i) = 0.0_dp
            else
               c%value(i) = a%value(mod(i-1,a%size())+1) / bj%value(1)
            end if
         end select
      end do
      if (present(status)) status = stat_i
   end function rate_binary

   function discount_multiply(a, b, status) result(c)
      type(discount_factor_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(discount_factor_t) :: c
      integer :: i, n, ia, ib
      n = max(a%size(), b%size())
      allocate(c%value(n), c%start_date(n), c%end_date(n))
      do i = 1, n
         ia = mod(i-1,a%size()) + 1
         ib = mod(i-1,b%size()) + 1
         if (a%end_date(ia) == b%start_date(ib)) then
            c%start_date(i) = a%start_date(ia)
            c%end_date(i) = b%end_date(ib)
         else if (b%end_date(ib) == a%start_date(ia)) then
            c%start_date(i) = b%start_date(ib)
            c%end_date(i) = a%end_date(ia)
         else
            c = empty_df()
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
         c%value(i) = a%value(ia) * b%value(ib)
      end do
      if (present(status)) status = FM_OK
   end function discount_multiply

   function discount_divide(numer, denom, status) result(c)
      type(discount_factor_t), intent(in) :: numer, denom
      integer, intent(out), optional :: status
      type(discount_factor_t) :: c
      integer :: i, n, ia, ib
      n = max(numer%size(), denom%size())
      allocate(c%value(n), c%start_date(n), c%end_date(n))
      do i = 1, n
         ia = mod(i-1,numer%size()) + 1
         ib = mod(i-1,denom%size()) + 1
         if (numer%start_date(ia) /= denom%start_date(ib) .or. &
             abs(denom%value(ib)) <= tiny(1.0_dp)) then
            c = empty_df()
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
         c%value(i) = numer%value(ia) / denom%value(ib)
         c%start_date(i) = denom%end_date(ib)
         c%end_date(i) = numer%end_date(ia)
         if (c%start_date(i) > c%end_date(i)) then
            c = empty_df()
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
      end do
      if (present(status)) status = FM_OK
   end function discount_divide

   pure integer function rate_size(self) result(value)
      class(interest_rate_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function rate_size

   pure integer function df_size(self) result(value)
      class(discount_factor_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function df_size

   subroutine validate_rate(rate, status)
      type(interest_rate_t), intent(in) :: rate
      integer, intent(out), optional :: status
      integer :: i
      if (rate%size() == 0 .or. size(rate%compounding) /= rate%size() .or. &
          size(rate%day_basis) /= rate%size()) then
         if (present(status)) status = FM_SIZE_MISMATCH
         return
      end if
      do i = 1, rate%size()
         if (.not. is_valid_compounding(rate%compounding(i)) .or. &
             .not. valid_day_basis(rate%day_basis(i))) then
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
      end do
      if (present(status)) status = FM_OK
   end subroutine validate_rate

   subroutine validate_df(df, status)
      type(discount_factor_t), intent(in) :: df
      integer, intent(out), optional :: status
      if (df%size() == 0 .or. size(df%start_date) /= df%size() .or. &
          size(df%end_date) /= df%size() .or. any(df%value <= 0.0_dp) .or. &
          any(df%start_date > df%end_date)) then
         if (present(status)) status = FM_INVALID_ARGUMENT
      else
         if (present(status)) status = FM_OK
      end if
   end subroutine validate_df

   pure logical function valid_day_basis(value) result(valid)
      character(len=*), intent(in) :: value
      character(len=12) :: b
      b = lower12(value)
      valid = trim(b) == 'act/365' .or. trim(b) == 'act/360' .or. &
         trim(b) == 'act/act' .or. trim(b) == '30e/360' .or. &
         trim(b) == '30/360us'
   end function valid_day_basis

   pure logical function is_continuous(compounding) result(value)
      real(dp), intent(in) :: compounding
      value = compounding >= 0.5_dp * huge(1.0_dp)
   end function is_continuous

   function empty_df() result(df)
      type(discount_factor_t) :: df
      allocate(df%value(0), df%start_date(0), df%end_date(0))
   end function empty_df

   pure function lower12(text) result(out)
      character(len=*), intent(in) :: text
      character(len=12) :: out
      integer :: i, k
      out = ''
      do i = 1, min(12, len_trim(text))
         k = iachar(text(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) then
            out(i:i) = achar(k + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower12

end module fmbasics_rates
