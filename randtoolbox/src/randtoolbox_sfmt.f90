! Scalar modern-Fortran translation of the SFMT implementation bundled with
! randtoolbox.  Upstream SFMT is BSD-3-Clause; see LICENSE.
module randtoolbox_sfmt
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox_base, only : u32_mod, u32, two32
   use randtoolbox_sfmt_params, only : get_sfmt_params, sfmt_supported
   implicit none
   private
   integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64)
   integer, save :: next_pset(6)=[1,1,1,1,1,1]

   type, public :: sfmt_rng
      integer :: mexp=19937, nblock=0, n32=0, idx=1
      integer :: pos1=0, sl1=0, sl2=0, sr1=0, sr2=0, pset=1
      integer(int64) :: masks(4)=0_int64, parity(4)=0_int64
      integer(int64), allocatable :: state(:)
      logical :: initialized=.false.
   contains
      procedure, public :: init => sfmt_init
      procedure, public :: next_uint32 => sfmt_next_uint32
      procedure, public :: next => sfmt_next
      procedure, public :: fill => sfmt_fill
      procedure, public :: fill_matrix => sfmt_fill_matrix
      procedure, public :: get_state => sfmt_get_state
   end type sfmt_rng

   public :: sfmt_random, reset_sfmt_parameter_sets, sfmt_supported
contains
   subroutine sfmt_init(this,mexp,seed,param_set)
      class(sfmt_rng), intent(inout) :: this
      integer, intent(in) :: mexp
      integer(int64), intent(in) :: seed
      integer, intent(in), optional :: param_set
      integer :: i, ps
      integer(int64) :: prev
      if(.not.sfmt_supported(mexp)) error stop 'randtoolbox: unsupported SFMT exponent'
      ps=1; if(present(param_set)) ps=param_set
      if(ps<1 .or. ps>32) error stop 'randtoolbox: SFMT parameter set must be 1..32'
      if(mexp>=44497) ps=1
      this%mexp=mexp; this%pset=ps
      this%nblock=mexp/128+1; this%n32=4*this%nblock
      call get_sfmt_params(mexp,ps,this%pos1,this%sl1,this%sl2,this%sr1,this%sr2,this%masks,this%parity)
      if(allocated(this%state)) deallocate(this%state)
      allocate(this%state(this%n32))
      this%state(1)=u32(seed)
      do i=2,this%n32
         prev=this%state(i-1)
         this%state(i)=u32(1812433253_int64*ieor(prev,shiftr(prev,30))+int(i-1,int64))
      end do
      call period_certification(this)
      this%idx=this%n32+1
      this%initialized=.true.
   end subroutine sfmt_init

   subroutine period_certification(this)
      class(sfmt_rng), intent(inout) :: this
      integer(int64) :: inner, work
      integer :: i,j,s
      inner=0_int64
      do i=1,4
         inner=ieor(inner,iand(this%state(i),this%parity(i)))
      end do
      s=16
      do while(s>0)
         inner=ieor(inner,shiftr(inner,s)); s=shiftr(s,1)
      end do
      if(iand(inner,1_int64)==1_int64) return
      do i=1,4
         work=1_int64
         do j=0,31
            if(iand(work,this%parity(i))/=0_int64) then
               this%state(i)=ieor(this%state(i),work)
               return
            end if
            work=shiftl(work,1)
         end do
      end do
   end subroutine period_certification

   function sfmt_next_uint32(this) result(r)
      class(sfmt_rng), intent(inout) :: this
      integer(int64) :: r
      if(.not.this%initialized) error stop 'randtoolbox: SFMT generator not initialized'
      if(this%idx>this%n32) then
         call gen_rand_all(this)
         this%idx=1
      end if
      r=this%state(this%idx)
      this%idx=this%idx+1
   end function sfmt_next_uint32

   real(real64) function sfmt_next(this) result(x)
      class(sfmt_rng), intent(inout) :: this
      x=(real(this%next_uint32(),real64)+0.5_real64)/two32
   end function sfmt_next

   subroutine sfmt_fill(this,x)
      class(sfmt_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:)
      integer :: i
      do i=1,size(x); x(i)=this%next(); end do
   end subroutine sfmt_fill

   subroutine sfmt_fill_matrix(this,x)
      class(sfmt_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:,:)
      integer :: i,j
      do j=1,size(x,2)
         do i=1,size(x,1)
            x(i,j)=this%next()
         end do
      end do
   end subroutine sfmt_fill_matrix

   subroutine sfmt_get_state(this,state,index)
      class(sfmt_rng), intent(in) :: this
      integer(int64), intent(out) :: state(:)
      integer, intent(out), optional :: index
      if(size(state)/=this%n32) error stop 'randtoolbox: SFMT state size mismatch'
      state=this%state
      if(present(index)) index=this%idx
   end subroutine sfmt_get_state

   subroutine gen_rand_all(this)
      class(sfmt_rng), intent(inout) :: this
      integer(int64) :: r1(4),r2(4),a(4),b(4),r(4)
      integer :: i, bi
      r1=get_block(this,this%nblock-1)
      r2=get_block(this,this%nblock)
      do i=1,this%nblock
         bi=i+this%pos1
         if(bi>this%nblock) bi=bi-this%nblock
         a=get_block(this,i); b=get_block(this,bi)
         call do_recursion(this,a,b,r1,r2,r)
         call put_block(this,i,r)
         r1=r2; r2=r
      end do
   end subroutine gen_rand_all

   pure function get_block(this,i) result(v)
      class(sfmt_rng), intent(in) :: this
      integer, intent(in) :: i
      integer(int64) :: v(4)
      v=this%state(4*(i-1)+1:4*i)
   end function get_block

   subroutine put_block(this,i,v)
      class(sfmt_rng), intent(inout) :: this
      integer, intent(in) :: i
      integer(int64), intent(in) :: v(4)
      this%state(4*(i-1)+1:4*i)=v
   end subroutine put_block

   subroutine do_recursion(this,a,b,c,d,r)
      class(sfmt_rng), intent(in) :: this
      integer(int64), intent(in) :: a(4),b(4),c(4),d(4)
      integer(int64), intent(out) :: r(4)
      integer(int64) :: x(4),y(4),t
      integer :: k
      call lshift128(a,this%sl2,x)
      call rshift128(c,this%sr2,y)
      do k=1,4
         t=ieor(a(k),x(k))
         t=ieor(t,iand(shiftr(b(k),this%sr1),this%masks(k)))
         t=ieor(t,y(k))
         t=ieor(t,shiftl(d(k),this%sl1))
         r(k)=iand(t,mask32)
      end do
   end subroutine do_recursion

   pure subroutine rshift128(v,shift,out)
      integer(int64), intent(in) :: v(4)
      integer, intent(in) :: shift
      integer(int64), intent(out) :: out(4)
      integer(int64) :: th,tl,oh,ol
      integer :: s
      s=8*shift
      th=ior(shiftl(v(4),32),v(3))
      tl=ior(shiftl(v(2),32),v(1))
      oh=shiftr(th,s)
      ol=shiftr(tl,s)
      if(s>0) ol=ior(ol,shiftl(th,64-s))
      out(1)=iand(ol,mask32); out(2)=iand(shiftr(ol,32),mask32)
      out(3)=iand(oh,mask32); out(4)=iand(shiftr(oh,32),mask32)
   end subroutine rshift128

   pure subroutine lshift128(v,shift,out)
      integer(int64), intent(in) :: v(4)
      integer, intent(in) :: shift
      integer(int64), intent(out) :: out(4)
      integer(int64) :: th,tl,oh,ol
      integer :: s
      s=8*shift
      th=ior(shiftl(v(4),32),v(3))
      tl=ior(shiftl(v(2),32),v(1))
      oh=shiftl(th,s)
      ol=shiftl(tl,s)
      if(s>0) oh=ior(oh,shiftr(tl,64-s))
      out(1)=iand(ol,mask32); out(2)=iand(shiftr(ol,32),mask32)
      out(3)=iand(oh,mask32); out(4)=iand(shiftr(oh,32),mask32)
   end subroutine lshift128

   function sfmt_random(n,dim,mexp,seed,use_parameter_sets) result(x)
      integer, intent(in) :: n,dim,mexp
      integer(int64), intent(in) :: seed
      logical, intent(in), optional :: use_parameter_sets
      real(real64), allocatable :: x(:,:)
      type(sfmt_rng) :: g
      logical :: ups
      integer :: slot,ps
      if(n<0 .or. dim<1) error stop 'randtoolbox: invalid SFMT dimensions'
      ups=.true.; if(present(use_parameter_sets)) ups=use_parameter_sets
      slot=mexp_slot(mexp)
      if(slot>0) then
         if(.not.ups) next_pset(slot)=1
         ps=next_pset(slot)
         next_pset(slot)=modulo(next_pset(slot),32)+1
      else
         ps=1
      end if
      call g%init(mexp,seed,ps)
      allocate(x(n,dim)); call g%fill_matrix(x)
   end function sfmt_random

   pure integer function mexp_slot(mexp) result(slot)
      integer,intent(in)::mexp
      select case(mexp)
      case(607);slot=1
      case(1279);slot=2
      case(2281);slot=3
      case(4253);slot=4
      case(11213);slot=5
      case(19937);slot=6
      case default;slot=0
      end select
   end function mexp_slot

   subroutine reset_sfmt_parameter_sets()
      next_pset=1
   end subroutine reset_sfmt_parameter_sets
end module randtoolbox_sfmt
