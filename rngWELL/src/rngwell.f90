! Modern Fortran translation of the computational code in rngWELL 0.10-10.
!
! rngWELL wrapper code: BSD 3-Clause; see LICENSE.
! WELL recurrence source notice: see LICENSES/WELL-SOURCE-NOTICE.txt.
module rngwell
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   use rngwell_bitops, only : lxor, rxor, lshift32, rshift32, mat2, mat4neg, mat5, &
      u32_to_real, u32_to_i64, i64_to_u32
   implicit none
   private

   integer, parameter, public :: max_well_state = 1391

   type, public :: well_rng
      private
      character(len=6) :: variant = ''
      integer(int32), allocatable :: state(:)
      integer :: idx = 1
   contains
      procedure, public :: init => well_init
      procedure, public :: seed => well_seed
      procedure, public :: next => well_next
      procedure, public :: next_uint32 => well_next_uint32
      procedure, public :: fill => well_fill
      procedure, public :: fill_matrix => well_fill_matrix
      procedure, public :: get_state => well_get_state
      procedure, public :: put_state => well_put_state
      procedure, public :: get_variant => well_get_variant
      procedure, public :: state_size => well_type_state_size
   end type well_rng

   public :: well_state_size, well_variant_supported, well_from_options, init_mt2002

contains

   pure integer function well_state_size(variant) result(n)
      character(len=*), intent(in) :: variant
      character(len=6) :: v
      v = canonical_variant(variant)
      select case (trim(v))
      case ('512a')
         n = 16
      case ('521a', '521b')
         n = 17
      case ('607a', '607b')
         n = 19
      case ('800a', '800b')
         n = 25
      case ('1024a', '1024b')
         n = 32
      case ('19937a', '19937b', '19937c')
         n = 624
      case ('21701a')
         n = 679
      case ('23209a', '23209b')
         n = 726
      case ('44497a', '44497b')
         n = 1391
      case default
         n = 0
      end select
   end function well_state_size

   pure logical function well_variant_supported(variant) result(ok)
      character(len=*), intent(in) :: variant
      ok = well_state_size(variant) > 0
   end function well_variant_supported

   pure character(len=6) function canonical_variant(variant) result(v)
      character(len=*), intent(in) :: variant
      integer :: i, c, j
      v = ''
      j = 0
      do i = 1, len_trim(variant)
         c = iachar(variant(i:i))
         if (variant(i:i) == ' ' .or. variant(i:i) == '-' .or. variant(i:i) == '_') cycle
         j = j + 1
         if (j > len(v)) exit
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            v(j:j) = achar(c + iachar('a') - iachar('A'))
         else
            v(j:j) = variant(i:i)
         end if
      end do
   end function canonical_variant

   subroutine well_init(this, variant, seed, state)
      class(well_rng), intent(inout) :: this
      character(len=*), intent(in) :: variant
      integer(int64), intent(in), optional :: seed
      integer(int32), intent(in), optional :: state(:)
      integer :: n
      integer(int64) :: clock_seed
      integer :: count

      this%variant = canonical_variant(variant)
      n = well_state_size(this%variant)
      if (n == 0) error stop 'rngWELL: unsupported WELL variant'
      if (allocated(this%state)) deallocate(this%state)
      allocate(this%state(n))
      this%idx = 1

      if (present(seed) .and. present(state)) then
         error stop 'rngWELL: provide seed or state, not both'
      else if (present(state)) then
         if (size(state) /= n) error stop 'rngWELL: wrong state length'
         this%state = state
      else if (present(seed)) then
         call this%seed(seed)
      else
         call system_clock(count=count)
         clock_seed = int(count, int64)
         call this%seed(clock_seed)
      end if
   end subroutine well_init

   subroutine well_seed(this, seed)
      class(well_rng), intent(inout) :: this
      integer(int64), intent(in) :: seed
      integer :: n

      n = well_state_size(this%variant)
      if (n == 0) error stop 'rngWELL: initialize variant before seeding'
      if (.not. allocated(this%state)) allocate(this%state(n))
      call init_mt2002(seed, this%state)
      this%idx = 1
   end subroutine well_seed

   subroutine init_mt2002(seed, state)
      integer(int64), intent(in) :: seed
      integer(int32), intent(out) :: state(:)
      integer :: i
      integer(int32) :: t
      integer(int64) :: u

      if (size(state) == 0) return
      state(1) = i64_to_u32(seed)
      do i = 2, size(state)
         t = ieor(state(i-1), shiftr(state(i-1), 30))
         u = 1812433253_int64 * u32_to_i64(t) + int(i-1, int64)
         state(i) = i64_to_u32(u)
      end do
   end subroutine init_mt2002

   integer function well_type_state_size(this) result(n)
      class(well_rng), intent(in) :: this
      n = well_state_size(this%variant)
   end function well_type_state_size

   function well_get_variant(this) result(v)
      class(well_rng), intent(in) :: this
      character(len=6) :: v
      v = this%variant
   end function well_get_variant

   subroutine well_get_state(this, state)
      class(well_rng), intent(in) :: this
      integer(int32), intent(out) :: state(:)
      integer :: n, n1

      if (.not. allocated(this%state)) error stop 'rngWELL: generator is not initialized'
      n = size(this%state)
      if (size(state) /= n) error stop 'rngWELL: wrong output state length'
      n1 = n - this%idx + 1
      state(1:n1) = this%state(this%idx:n)
      if (this%idx > 1) state(n1+1:n) = this%state(1:this%idx-1)
   end subroutine well_get_state

   subroutine well_put_state(this, variant, state)
      class(well_rng), intent(inout) :: this
      character(len=*), intent(in) :: variant
      integer(int32), intent(in) :: state(:)
      call this%init(variant, state=state)
   end subroutine well_put_state

   subroutine well_fill(this, x)
      class(well_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:)
      integer :: i
      do i = 1, size(x)
         x(i) = this%next()
      end do
   end subroutine well_fill

   subroutine well_fill_matrix(this, x)
      class(well_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:,:)
      integer :: i, j
      ! Matches rngWELL's column-major R matrix fill order.
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            x(i,j) = this%next()
         end do
      end do
   end subroutine well_fill_matrix

   function well_from_options(order, version, temper, seed) result(rng)
      integer, intent(in) :: order
      character(len=*), intent(in) :: version
      logical, intent(in), optional :: temper
      integer(int64), intent(in), optional :: seed
      type(well_rng) :: rng
      character(len=1) :: ver
      character(len=6) :: variant
      logical :: t
      integer(int64) :: s

      t = .false.
      if (present(temper)) t = temper
      ver = canonical_version(version)
      if (ver == '?') error stop 'rngWELL: version must be a or b'

      select case (order)
      case (512)
         if (ver /= 'a') error stop 'rngWELL: WELL512 has only version a'
         variant = '512a'
         if (t) error stop 'rngWELL: WELL512 cannot be tempered'
      case (521)
         variant = '521' // ver
         if (t) error stop 'rngWELL: WELL521 cannot be tempered'
      case (607)
         variant = '607' // ver
         if (t) error stop 'rngWELL: WELL607 cannot be tempered'
      case (800)
         variant = '800' // ver
         ! Upstream WELL2test accepts temper=.true. here but ignores it.
      case (1024)
         variant = '1024' // ver
         if (t) error stop 'rngWELL: WELL1024 cannot be tempered'
      case (19937)
         if (t) then
            variant = '19937c'
         else
            variant = '19937' // ver
         end if
      case (21701)
         if (ver /= 'a') error stop 'rngWELL: WELL21701 has only version a'
         variant = '21701a'
         ! Upstream accepts temper but ignores it.
      case (23209)
         variant = '23209' // ver
         ! Upstream accepts temper but ignores it.
      case (44497)
         if (t) then
            variant = '44497b'
         else
            ! WELL2test uses 44497a when temper is false, irrespective of version.
            variant = '44497a'
         end if
      case default
         error stop 'rngWELL: unsupported WELL order'
      end select

      if (present(seed)) then
         s = seed
         call rng%init(trim(variant), seed=s)
      else
         call rng%init(trim(variant))
      end if
   end function well_from_options

   pure character(len=1) function canonical_version(version) result(v)
      character(len=*), intent(in) :: version
      integer :: c
      if (len_trim(version) == 0) then
         v = '?'
         return
      end if
      c = iachar(version(1:1))
      if (c >= iachar('A') .and. c <= iachar('Z')) c = c + iachar('a') - iachar('A')
      select case (achar(c))
      case ('a', 'b')
         v = achar(c)
      case default
         v = '?'
      end select
   end function canonical_version

   real(real64) function well_next(this) result(u)
      class(well_rng), intent(inout) :: this
      integer(int32) :: y

      y = this%next_uint32()
      u = u32_to_real(y)
   end function well_next

   integer(int32) function well_next_uint32(this) result(y)
      class(well_rng), intent(inout) :: this

      if (.not. allocated(this%state)) error stop 'rngWELL: generator is not initialized'

      select case (trim(this%variant))
      case ('512a')
         y = step_512a(this%state, this%idx)
      case ('521a')
         y = step_521a(this%state, this%idx)
      case ('521b')
         y = step_521b(this%state, this%idx)
      case ('607a')
         y = step_607a(this%state, this%idx)
      case ('607b')
         y = step_607b(this%state, this%idx)
      case ('800a')
         y = step_800a(this%state, this%idx)
      case ('800b')
         y = step_800b(this%state, this%idx)
      case ('1024a')
         y = step_1024a(this%state, this%idx)
      case ('1024b')
         y = step_1024b(this%state, this%idx)
      case ('19937a')
         y = step_19937a(this%state, this%idx, .false.)
      case ('19937b')
         y = step_19937b(this%state, this%idx)
      case ('19937c')
         y = step_19937a(this%state, this%idx, .true.)
      case ('21701a')
         y = step_21701a(this%state, this%idx)
      case ('23209a')
         y = step_23209a(this%state, this%idx)
      case ('23209b')
         y = step_23209b(this%state, this%idx)
      case ('44497a')
         y = step_44497a(this%state, this%idx, .false.)
      case ('44497b')
         y = step_44497a(this%state, this%idx, .true.)
      case default
         error stop 'rngWELL: unsupported WELL variant'
      end select
   end function well_next_uint32

   pure integer function at(idx, offset, n) result(j)
      integer, intent(in) :: idx, offset, n
      j = modulo(idx - 1 + offset, n) + 1
   end function at

   pure integer(int32) function z0_masked(state, idx, p) result(z0)
      integer(int32), intent(in) :: state(:)
      integer, intent(in) :: idx, p
      integer(int32) :: masku, maskl
      integer :: n
      n = size(state)
      if (p == 0) then
         z0 = state(at(idx, -1, n))
      else
         masku = shiftr(not(0_int32), 32-p)
         maskl = not(masku)
         z0 = ior(iand(state(at(idx, -1, n)), maskl), iand(state(at(idx, -2, n)), masku))
      end if
   end function z0_masked

   integer(int32) function finish_step(state, idx, newv1, newv0) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32), intent(in) :: newv1, newv0
      integer :: n, prev
      n = size(state)
      prev = at(idx, -1, n)
      state(idx) = newv1
      state(prev) = newv0
      idx = prev
      y = state(idx)
   end function finish_step

   integer(int32) function step_512a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0, z1, z2, nv1, nv0
      integer :: n
      n = 16
      z0 = state(at(idx,-1,n))
      z1 = ieor(lxor(state(idx),16), lxor(state(at(idx,13,n)),15))
      z2 = rxor(state(at(idx,9,n)),11)
      nv1 = ieor(z1,z2)
      nv0 = ieor(ieor(lxor(z0,2), lxor(z1,18)), ieor(lshift32(z2,28), &
         mat4neg(5,int(z'da442d24',int32),nv1)))
      y = finish_step(state,idx,nv1,nv0)
   end function step_512a

   integer(int32) function step_521a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=17
      z0=z0_masked(state,idx,23)
      z1=ieor(lxor(state(idx),13),lxor(state(at(idx,13,n)),15))
      z2=ieor(state(at(idx,11,n)),lshift32(state(at(idx,10,n)),21))
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,13),ieor(rshift32(z1,1),rxor(nv1,11)))
      y=finish_step(state,idx,nv1,nv0)
   end function step_521a

   integer(int32) function step_521b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=17
      z0=z0_masked(state,idx,23)
      z1=ieor(lxor(state(idx),21),rxor(state(at(idx,11,n)),6))
      z2=lxor(state(at(idx,7,n)),13)
      nv1=ieor(z1,z2)
      nv0=ieor(rxor(z0,13),ieor(lshift32(z1,10),ieor(lshift32(z2,5),rxor(nv1,13))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_521b

   integer(int32) function step_607a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=19
      z0=z0_masked(state,idx,1)
      z1=ieor(rxor(state(idx),19),rxor(state(at(idx,16,n)),11))
      z2=ieor(lxor(state(at(idx,15,n)),14),state(at(idx,14,n)))
      nv1=ieor(z1,z2)
      nv0=ieor(rxor(z0,18),ieor(z1,lxor(nv1,5)))
      y=finish_step(state,idx,nv1,nv0)
   end function step_607a

   integer(int32) function step_607b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=19
      z0=z0_masked(state,idx,1)
      z1=ieor(lxor(state(idx),18),lxor(state(at(idx,16,n)),14))
      z2=rxor(state(at(idx,13,n)),18)
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,24),ieor(rxor(z1,5),lxor(z2,1)))
      y=finish_step(state,idx,nv1,nv0)
   end function step_607b

   integer(int32) function step_800a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=25
      z0=state(at(idx,-1,n))
      z1=ieor(state(idx),lxor(state(at(idx,14,n)),15))
      z2=ieor(rxor(state(at(idx,18,n)),10),lxor(state(at(idx,17,n)),11))
      nv1=ieor(z1,z2)
      nv0=ieor(rxor(z0,16),ieor(rshift32(z1,20),ieor(z2,lxor(nv1,28))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_800a

   integer(int32) function step_800b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=25
      z0=state(at(idx,-1,n))
      z1=ieor(lxor(state(idx),29),lshift32(state(at(idx,9,n)),14))
      z2=ieor(state(at(idx,4,n)),rshift32(state(at(idx,22,n)),19))
      nv1=ieor(z1,z2)
      nv0=ieor(z0,ieor(rxor(z1,10),ieor(mat2(int(z'd3e43ffd',int32),z2),lxor(nv1,25))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_800b

   integer(int32) function step_1024a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=32
      z0=state(at(idx,-1,n))
      z1=ieor(state(idx),rxor(state(at(idx,3,n)),8))
      z2=ieor(lxor(state(at(idx,24,n)),19),lxor(state(at(idx,10,n)),14))
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,11),ieor(lxor(z1,7),lxor(z2,13)))
      y=finish_step(state,idx,nv1,nv0)
   end function step_1024a

   integer(int32) function step_1024b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=32
      z0=state(at(idx,-1,n))
      z1=ieor(lxor(state(idx),21),rxor(state(at(idx,22,n)),17))
      z2=ieor(mat2(int(z'8bdcb91e',int32),state(at(idx,25,n))),rxor(state(at(idx,26,n)),15))
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,14),ieor(lxor(z1,21),z2))
      y=finish_step(state,idx,nv1,nv0)
   end function step_1024b

   integer(int32) function step_19937a(state, idx, tempered) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      logical, intent(in) :: tempered
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=624
      z0=z0_masked(state,idx,31)
      z1=ieor(lxor(state(idx),25),rxor(state(at(idx,70,n)),27))
      z2=ieor(rshift32(state(at(idx,179,n)),9),rxor(state(at(idx,449,n)),1))
      nv1=ieor(z1,z2)
      nv0=ieor(z0,ieor(lxor(z1,9),ieor(lxor(z2,21),rxor(nv1,21))))
      y=finish_step(state,idx,nv1,nv0)
      if (tempered) then
         y=ieor(y,iand(shiftl(y,7),int(z'e46e1700',int32)))
         y=ieor(y,iand(shiftl(y,15),int(z'9b868000',int32)))
      end if
   end function step_19937a

   integer(int32) function step_19937b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=624
      z0=z0_masked(state,idx,31)
      z1=ieor(rxor(state(idx),7),state(at(idx,203,n)))
      z2=ieor(rxor(state(at(idx,613,n)),12),lxor(state(at(idx,123,n)),10))
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,19),ieor(lshift32(z1,11),ieor(rxor(z2,4),lxor(nv1,10))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_19937b

   integer(int32) function step_21701a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=679
      z0=z0_masked(state,idx,27)
      z1=ieor(state(idx),lxor(state(at(idx,151,n)),26))
      z2=rxor(state(at(idx,327,n)),19)
      nv1=ieor(z1,z2)
      nv0=ieor(rxor(z0,27),ieor(lxor(z1,11),ieor( &
         mat5(15,int(z'86a9d87e',int32),int(z'ffffffef',int32),int(z'00200000',int32),z2), &
         lxor(nv1,16))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_21701a

   integer(int32) function step_23209a(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=726
      z0=z0_masked(state,idx,23)
      z1=ieor(rxor(state(idx),28),state(at(idx,667,n)))
      z2=ieor(rxor(state(at(idx,43,n)),18),rxor(state(at(idx,462,n)),3))
      nv1=ieor(z1,z2)
      nv0=ieor(rxor(z0,21),ieor(lxor(z1,17),ieor(lxor(z2,28),lxor(nv1,1))))
      y=finish_step(state,idx,nv1,nv0)
   end function step_23209a

   integer(int32) function step_23209b(state, idx) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=726
      z0=z0_masked(state,idx,23)
      z1=ieor(mat2(int(z'a8c296d1',int32),state(idx)),state(at(idx,610,n)))
      z2=ieor(mat5(15,int(z'5d6b45cc',int32),int(z'fffeffff',int32),int(z'00000002',int32), &
         state(at(idx,175,n))),lxor(state(at(idx,662,n)),24))
      nv1=ieor(z1,z2)
      nv0=ieor(lxor(z0,26),ieor(z1,rxor(nv1,16)))
      y=finish_step(state,idx,nv1,nv0)
   end function step_23209b

   integer(int32) function step_44497a(state, idx, tempered) result(y)
      integer(int32), intent(inout) :: state(:)
      integer, intent(inout) :: idx
      logical, intent(in) :: tempered
      integer(int32) :: z0,z1,z2,nv1,nv0
      integer :: n
      n=1391
      z0=z0_masked(state,idx,15)
      z1=ieor(lxor(state(idx),24),rxor(state(at(idx,23,n)),30))
      z2=ieor(lxor(state(at(idx,481,n)),10),lshift32(state(at(idx,229,n)),26))
      nv1=ieor(z1,z2)
      nv0=ieor(z0,ieor(rxor(z1,20),ieor( &
         mat5(9,int(z'b729fcec',int32),int(z'fbffffff',int32),int(z'00020000',int32),z2),nv1)))
      y=finish_step(state,idx,nv1,nv0)
      if (tempered) then
         y=ieor(y,iand(shiftl(y,7),int(z'93dd1400',int32)))
         y=ieor(y,iand(shiftl(y,15),int(z'fa118000',int32)))
      end if
   end function step_44497a

end module rngwell
