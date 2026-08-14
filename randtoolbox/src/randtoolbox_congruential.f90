module randtoolbox_congruential
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox_base, only : two64, u64_real, mul_add_low64, system_seed
   implicit none
   private

   type, public :: congru_rng
      integer(int64) :: modulus = 2147483647_int64
      integer(int64) :: multiplier = 16807_int64
      integer(int64) :: increment = 0_int64
      integer(int64) :: state = 1_int64
      integer(int64) :: mask = 0_int64
      logical :: power2_mode = .false.
      logical :: mod2_64_mode = .false.
   contains
      procedure, public :: init => congru_init
      procedure, public :: next => congru_next
      procedure, public :: fill => congru_fill
      procedure, public :: fill_matrix => congru_fill_matrix
      procedure, public :: get_state => congru_get_state
   end type congru_rng

   public :: modular_multiply

contains
   subroutine congru_init(this, seed, modulus, multiplier, increment)
      class(congru_rng), intent(inout) :: this
      integer(int64), intent(in), optional :: seed, modulus, multiplier, increment
      integer(int64) :: s, m

      if (present(modulus)) this%modulus = modulus
      if (present(multiplier)) this%multiplier = multiplier
      if (present(increment)) this%increment = increment
      s = system_seed(); if (present(seed)) s = seed

      if (this%modulus < 0_int64) error stop 'randtoolbox: modulus must be positive or zero for 2^64'
      if (this%multiplier <= 0_int64) error stop 'randtoolbox: multiplier must be positive'
      if (this%increment < 0_int64) error stop 'randtoolbox: increment must be nonnegative'

      if (this%modulus == 0_int64) then
         this%mod2_64_mode = .true.
         this%power2_mode = .false.
         this%mask = -1_int64
         this%state = s
      else
         this%mod2_64_mode = .false.
         if (this%multiplier >= this%modulus) error stop 'randtoolbox: multiplier must be less than modulus'
         if (this%increment >= this%modulus) error stop 'randtoolbox: increment must be less than modulus'
         m = this%modulus
         this%power2_mode = (iand(m, m-1_int64) == 0_int64)
         if (this%power2_mode) then
            this%mask = m - 1_int64
         else
            this%mask = 0_int64
         end if
         this%state = modulo(s, this%modulus)
      end if
   end subroutine congru_init

   real(real64) function congru_next(this) result(x)
      class(congru_rng), intent(inout) :: this
      integer(int64) :: t
      if (this%mod2_64_mode) then
         this%state = mul_add_low64(this%multiplier, this%state, this%increment)
         x = u64_real(this%state) / two64
      else if (this%power2_mode) then
         t = mul_add_low64(this%multiplier, this%state, this%increment)
         this%state = iand(t, this%mask)
         x = real(this%state, real64) / real(this%modulus, real64)
      else
         this%state = add_mod(modular_multiply(this%multiplier, this%state, this%modulus), &
                              this%increment, this%modulus)
         x = real(this%state, real64) / real(this%modulus, real64)
      end if
      ! Match randtoolbox's standalone congruRand convention.
      if (this%state == 0_int64) x = 1.0_real64
   end function congru_next

   pure integer(int64) function modular_multiply(a,b,m) result(r)
      integer(int64), intent(in) :: a,b,m
      integer(int64) :: x, y
      if (m <= 0_int64) error stop 'randtoolbox: modular_multiply requires positive modulus'
      x = modulo(a,m); y = modulo(b,m); r = 0_int64
      do while (y > 0_int64)
         if (btest(y,0)) r = add_mod(r,x,m)
         y = shiftr(y,1)
         if (y > 0_int64) x = add_mod(x,x,m)
      end do
   end function modular_multiply

   pure integer(int64) function add_mod(a,b,m) result(r)
      integer(int64), intent(in) :: a,b,m
      ! a,b in [0,m); avoid signed overflow.
      if (a >= m-b) then
         r = a - (m-b)
      else
         r = a+b
      end if
   end function add_mod

   subroutine congru_fill(this, x)
      class(congru_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:)
      integer :: i
      do i=1,size(x); x(i)=this%next(); end do
   end subroutine congru_fill

   subroutine congru_fill_matrix(this, x)
      class(congru_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:,:)
      integer :: i,j
      ! Upstream congruRand loops rows outside columns: each row receives
      ! consecutive values across dimensions.
      do i=1,size(x,1)
         do j=1,size(x,2)
            x(i,j)=this%next()
         end do
      end do
   end subroutine congru_fill_matrix

   pure integer(int64) function congru_get_state(this) result(s)
      class(congru_rng), intent(in) :: this
      s=this%state
   end function congru_get_state
end module randtoolbox_congruential
