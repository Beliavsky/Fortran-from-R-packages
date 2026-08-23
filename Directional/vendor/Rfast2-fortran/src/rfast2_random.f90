module rfast2_random
   use, intrinsic :: iso_fortran_env, only : int64
   use rfast_special, only : dp, pi
   implicit none
   private

   integer(int64), parameter :: mask16 = int(z'FFFF', int64)
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)
   integer(int64), parameter :: base16 = 65536_int64
   integer(int64), parameter :: mul_limb(0:3) = [32557_int64, 19605_int64, 62509_int64, 22609_int64]

   type :: pcg32_state
      integer(int64) :: state = 0_int64
      integer(int64) :: inc = 1442695040888963407_int64
   contains
      procedure :: seed => pcg_seed
      procedure :: next_u32 => pcg_next_u32
      procedure :: uniform_open => pcg_uniform_open
      procedure :: uniform_closed => pcg_uniform_closed
   end type pcg32_state

   type(pcg32_state), save :: rng
   logical, save :: initialized = .false.
   logical, save :: has_spare = .false.
   real(dp), save :: spare = 0.0_dp

   public :: set_seed, runif, sample_int, sample_real
   public :: rbeta_fast, rexp_fast, rchisq_fast, rgamma_fast, rgeom_fast, rcauchy_fast, rt_fast
   public :: rbeta1

contains

   pure function pcg_muladd(x, inc) result(y)
      integer(int64), intent(in) :: x, inc
      integer(int64) :: y
      integer(int64) :: s(0:3), q(0:3), r(0:3), total, carry
      integer :: k, j

      do k = 0, 3
         s(k) = iand(shiftr(x, 16*k), mask16)
         q(k) = iand(shiftr(inc, 16*k), mask16)
      end do
      carry = 0_int64
      do k = 0, 3
         total = carry + q(k)
         do j = 0, k
            total = total + s(j) * mul_limb(k-j)
         end do
         r(k) = modulo(total, base16)
         carry = total / base16
      end do
      y = r(0)
      y = ior(y, shiftl(r(1),16))
      y = ior(y, shiftl(r(2),32))
      y = ior(y, shiftl(r(3),48))
   end function pcg_muladd

   subroutine pcg_seed(self, seed_value, stream)
      class(pcg32_state), intent(inout) :: self
      integer(int64), intent(in) :: seed_value
      integer(int64), intent(in), optional :: stream

      self%state = seed_value
      if (present(stream)) then
         self%inc = ior(stream,1_int64)
      else
         self%inc = ior(seed_value,1_int64)
      end if
   end subroutine pcg_seed

   function pcg_next_u32(self) result(out)
      class(pcg32_state), intent(inout) :: self
      integer(int64) :: out, oldstate, xorshifted
      integer :: rot

      oldstate = self%state
      self%state = pcg_muladd(oldstate,self%inc)
      xorshifted = iand(shiftr(ieor(shiftr(oldstate,18),oldstate),27),mask32)
      rot = int(iand(shiftr(oldstate,59),31_int64))
      if (rot == 0) then
         out = xorshifted
      else
         out = iand(ior(shiftr(xorshifted,rot),shiftl(xorshifted,32-rot)),mask32)
      end if
   end function pcg_next_u32

   function pcg_uniform_open(self) result(u)
      class(pcg32_state), intent(inout) :: self
      real(dp) :: u

      u = (real(self%next_u32(),dp) + 0.5_dp) / 4294967296.0_dp
   end function pcg_uniform_open

   function pcg_uniform_closed(self) result(u)
      class(pcg32_state), intent(inout) :: self
      real(dp) :: u

      u = real(self%next_u32(),dp) / 4294967295.0_dp
   end function pcg_uniform_closed

   subroutine ensure_rng()
      integer :: count
      if (initialized) return
      call system_clock(count)
      call set_seed(int(count,int64))
   end subroutine ensure_rng

   subroutine set_seed(seed_value)
      integer(int64), intent(in) :: seed_value

      call rng%seed(seed_value,seed_value)
      has_spare = .false.
      initialized = .true.
   end subroutine set_seed

   real(dp) function normal_standard() result(z)
      real(dp) :: u1, u2, r

      call ensure_rng()
      if (has_spare) then
         z = spare
         has_spare = .false.
         return
      end if
      u1 = rng%uniform_open()
      u2 = rng%uniform_open()
      r = sqrt(-2.0_dp*log(u1))
      z = r*cos(2.0_dp*pi*u2)
      spare = r*sin(2.0_dp*pi*u2)
      has_spare = .true.
   end function normal_standard

   recursive real(dp) function gamma_one(shape, rate) result(x)
      real(dp), intent(in) :: shape, rate
      real(dp) :: d, c, z, v, u

      if (shape <= 0.0_dp .or. rate <= 0.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      if (shape < 1.0_dp) then
         x = gamma_one(shape+1.0_dp,rate) * rng%uniform_open()**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         z = normal_standard()
         v = 1.0_dp + c*z
         if (v <= 0.0_dp) cycle
         v = v*v*v
         u = rng%uniform_open()
         if (u < 1.0_dp - 0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
      end do
      x = d*v/rate
   end function gamma_one

   function runif(n, minval, maxval) result(x)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: minval, maxval
      real(dp), allocatable :: x(:)
      real(dp) :: a, b
      integer :: i

      call ensure_rng()
      a = 0.0_dp
      b = 1.0_dp
      if (present(minval)) a = minval
      if (present(maxval)) b = maxval
      allocate(x(max(0,n)))
      do i = 1, size(x)
         x(i) = a + (b-a)*rng%uniform_closed()
      end do
   end function runif

   function sample_int(n, size_out, replace) result(x)
      integer, intent(in) :: n, size_out
      logical, intent(in), optional :: replace
      integer, allocatable :: x(:), pool(:)
      logical :: repl
      integer :: i, j, m, tmp

      call ensure_rng()
      repl = .false.
      if (present(replace)) repl = replace
      if (n <= 0 .or. size_out < 0 .or. (.not. repl .and. size_out > n)) then
         allocate(x(0))
         return
      end if
      allocate(x(size_out))
      if (repl) then
         do i = 1, size_out
            x(i) = 1 + int(rng%uniform_open()*real(n,dp))
            x(i) = min(n,x(i))
         end do
      else
         allocate(pool(n))
         pool = [(i,i=1,n)]
         m = n
         do i = 1, size_out
            j = 1 + int(rng%uniform_open()*real(m,dp))
            j = min(m,j)
            x(i) = pool(j)
            tmp = pool(j)
            pool(j) = pool(m)
            pool(m) = tmp
            m = m - 1
         end do
      end if
   end function sample_int

   function sample_real(values, size_out, replace) result(x)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: size_out
      logical, intent(in), optional :: replace
      real(dp), allocatable :: x(:)
      integer, allocatable :: idx(:)

      if (present(replace)) then
         idx = sample_int(size(values),size_out,replace)
      else
         idx = sample_int(size(values),size_out)
      end if
      allocate(x(size(idx)))
      if (size(idx) > 0) x = values(idx)
   end function sample_real

   function rgamma_fast(n, shape, rate) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: shape
      real(dp), intent(in), optional :: rate
      real(dp), allocatable :: x(:)
      real(dp) :: r
      integer :: i

      call ensure_rng()
      r = 1.0_dp
      if (present(rate)) r = rate
      allocate(x(max(0,n)))
      do i = 1, size(x)
         x(i) = gamma_one(shape,r)
      end do
   end function rgamma_fast

   function rbeta_fast(n, alpha, beta) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: alpha, beta
      real(dp), allocatable :: x(:)
      real(dp) :: a, b
      integer :: i

      call ensure_rng()
      allocate(x(max(0,n)))
      do i = 1, size(x)
         a = gamma_one(alpha,1.0_dp)
         b = gamma_one(beta,1.0_dp)
         x(i) = a/(a+b)
      end do
   end function rbeta_fast

   function rbeta1(n, alpha) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), allocatable :: x(:)
      real(dp) :: a, b
      integer :: i

      call ensure_rng()
      allocate(x(max(0,n)))
      do i = 1, size(x)
         a = gamma_one(alpha,1.0_dp)
         b = gamma_one(alpha,1.0_dp)
         x(i) = a/(a+b)
      end do
   end function rbeta1

   function rexp_fast(n, rate) result(x)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: rate
      real(dp), allocatable :: x(:)
      real(dp) :: r
      integer :: i

      call ensure_rng()
      r = 1.0_dp
      if (present(rate)) r = rate
      allocate(x(max(0,n)))
      do i = 1, size(x)
         x(i) = -log(rng%uniform_open())/r
      end do
   end function rexp_fast

   function rchisq_fast(n, df) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: df
      real(dp), allocatable :: x(:)
      integer :: i

      call ensure_rng()
      allocate(x(max(0,n)))
      do i = 1, size(x)
         x(i) = gamma_one(0.5_dp*df,0.5_dp)
      end do
   end function rchisq_fast

   function rgeom_fast(n, prob) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: prob
      integer, allocatable :: x(:)
      integer :: i

      call ensure_rng()
      allocate(x(max(0,n)))
      do i = 1, size(x)
         if (prob >= 1.0_dp) then
            x(i) = 0
         else
            x(i) = floor(log(rng%uniform_open())/log(1.0_dp-prob))
         end if
      end do
   end function rgeom_fast

   function rcauchy_fast(n, location, scale) result(x)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: location, scale
      real(dp), allocatable :: x(:)
      real(dp) :: loc, sc
      integer :: i

      call ensure_rng()
      loc = 0.0_dp
      sc = 1.0_dp
      if (present(location)) loc = location
      if (present(scale)) sc = scale
      allocate(x(max(0,n)))
      do i = 1, size(x)
         x(i) = loc + sc*tan(pi*(rng%uniform_open()-0.5_dp))
      end do
   end function rcauchy_fast

   function rt_fast(n, df, ncp) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: df
      real(dp), intent(in), optional :: ncp
      real(dp), allocatable :: x(:)
      real(dp) :: shift, cs
      integer :: i

      call ensure_rng()
      shift = 0.0_dp
      if (present(ncp)) shift = ncp
      allocate(x(max(0,n)))
      do i = 1, size(x)
         cs = gamma_one(0.5_dp*df,0.5_dp)
         x(i) = (normal_standard()+shift)/sqrt(cs/df)
      end do
   end function rt_fast

end module rfast2_random
