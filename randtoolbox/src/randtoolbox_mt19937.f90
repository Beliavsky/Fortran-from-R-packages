! Modern Fortran translation of mt19937ar.c bundled with randtoolbox.
! Copyright and BSD-3-Clause terms: see LICENSES/MT19937-BSD-3-Clause.txt.
module randtoolbox_mt19937
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox_base, only : u32, two32
   implicit none
   private
   integer, parameter :: nstate=624, mlag=397
   integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64)
   integer(int64), parameter :: matrix_a=int(z'9908B0DF',int64)
   integer(int64), parameter :: upper_mask=int(z'80000000',int64)
   integer(int64), parameter :: lower_mask=int(z'7FFFFFFF',int64)

   type, public :: mt19937_rng
      integer(int64) :: state(nstate)=0_int64
      integer :: idx=nstate+2
   contains
      procedure, public :: init => mt_init_seed
      procedure, public :: init_array => mt_init_array
      procedure, public :: next_uint32 => mt_next_uint32
      procedure, public :: next32 => mt_next32
      procedure, public :: next53 => mt_next53
      procedure, public :: fill => mt_fill
      procedure, public :: fill_matrix => mt_fill_matrix
      procedure, public :: get_state => mt_get_state
      procedure, public :: put_state => mt_put_state
   end type mt19937_rng
contains
   subroutine mt_init_seed(this,seed)
      class(mt19937_rng), intent(inout) :: this
      integer(int64), intent(in) :: seed
      integer :: i
      integer(int64) :: p
      this%state(1)=u32(seed)
      do i=2,nstate
         p=this%state(i-1)
         this%state(i)=u32(1812433253_int64*ieor(p,shiftr(p,30))+int(i-1,int64))
      end do
      this%idx=nstate+1
   end subroutine mt_init_seed

   subroutine mt_init_array(this,key)
      class(mt19937_rng), intent(inout) :: this
      integer(int64), intent(in) :: key(:)
      integer :: i,j,k
      integer(int64) :: p,t
      if(size(key)<1) error stop 'randtoolbox: MT19937 key array must be nonempty'
      call this%init(19650218_int64)
      i=2; j=1
      do k=1,max(nstate,size(key))
         p=this%state(i-1)
         t=u32(ieor(p,shiftr(p,30))*1664525_int64)
         this%state(i)=u32(ieor(this%state(i),t)+u32(key(j))+int(j-1,int64))
         i=i+1; j=j+1
         if(i>nstate) then
            this%state(1)=this%state(nstate); i=2
         end if
         if(j>size(key)) j=1
      end do
      do k=1,nstate-1
         p=this%state(i-1)
         t=u32(ieor(p,shiftr(p,30))*1566083941_int64)
         this%state(i)=u32(ieor(this%state(i),t)-int(i-1,int64))
         i=i+1
         if(i>nstate) then
            this%state(1)=this%state(nstate); i=2
         end if
      end do
      this%state(1)=upper_mask
      this%idx=nstate+1
   end subroutine mt_init_array

   function mt_next_uint32(this) result(y)
      class(mt19937_rng), intent(inout) :: this
      integer(int64) :: y
      integer :: k
      if(this%idx>nstate) then
         if(this%idx==nstate+2) call this%init(5489_int64)
         do k=1,nstate-mlag
            y=ior(iand(this%state(k),upper_mask),iand(this%state(k+1),lower_mask))
            this%state(k)=u32(ieor(this%state(k+mlag),shiftr(y,1)))
            if(iand(y,1_int64)/=0_int64) this%state(k)=ieor(this%state(k),matrix_a)
         end do
         do k=nstate-mlag+1,nstate-1
            y=ior(iand(this%state(k),upper_mask),iand(this%state(k+1),lower_mask))
            this%state(k)=u32(ieor(this%state(k+mlag-nstate),shiftr(y,1)))
            if(iand(y,1_int64)/=0_int64) this%state(k)=ieor(this%state(k),matrix_a)
         end do
         y=ior(iand(this%state(nstate),upper_mask),iand(this%state(1),lower_mask))
         this%state(nstate)=u32(ieor(this%state(mlag),shiftr(y,1)))
         if(iand(y,1_int64)/=0_int64) this%state(nstate)=ieor(this%state(nstate),matrix_a)
         this%idx=1
      end if
      y=this%state(this%idx); this%idx=this%idx+1
      y=ieor(y,shiftr(y,11))
      y=u32(ieor(y,iand(shiftl(y,7),int(z'9D2C5680',int64))))
      y=u32(ieor(y,iand(shiftl(y,15),int(z'EFC60000',int64))))
      y=u32(ieor(y,shiftr(y,18)))
   end function mt_next_uint32

   real(real64) function mt_next32(this) result(x)
      class(mt19937_rng), intent(inout) :: this
      x=(real(this%next_uint32(),real64)+0.5_real64)/two32
   end function mt_next32

   real(real64) function mt_next53(this) result(x)
      class(mt19937_rng), intent(inout) :: this
      integer(int64) :: a,b
      a=shiftr(this%next_uint32(),5); b=shiftr(this%next_uint32(),6)
      x=(real(a,real64)*67108864.0_real64+real(b,real64))/9007199254740992.0_real64
   end function mt_next53

   subroutine mt_fill(this,x,resolution)
      class(mt19937_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:)
      integer, intent(in), optional :: resolution
      integer :: i,r
      r=32; if(present(resolution)) r=resolution
      if(r/=32 .and. r/=53) error stop 'randtoolbox: MT resolution must be 32 or 53'
      do i=1,size(x)
         if(r==32) then; x(i)=this%next32(); else; x(i)=this%next53(); end if
      end do
   end subroutine mt_fill

   subroutine mt_fill_matrix(this,x,resolution)
      class(mt19937_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:,:)
      integer, intent(in), optional :: resolution
      integer :: i,j,r
      r=32; if(present(resolution)) r=resolution
      if(r/=32 .and. r/=53) error stop 'randtoolbox: MT resolution must be 32 or 53'
      do j=1,size(x,2); do i=1,size(x,1)
         if(r==32) then; x(i,j)=this%next32(); else; x(i,j)=this%next53(); end if
      end do; end do
   end subroutine mt_fill_matrix

   subroutine mt_get_state(this,state,index)
      class(mt19937_rng), intent(in) :: this
      integer(int64), intent(out) :: state(:)
      integer, intent(out), optional :: index
      if(size(state)/=nstate) error stop 'randtoolbox: MT state size mismatch'
      state=this%state; if(present(index)) index=this%idx
   end subroutine mt_get_state

   subroutine mt_put_state(this,state,index)
      class(mt19937_rng), intent(inout) :: this
      integer(int64), intent(in) :: state(:)
      integer, intent(in) :: index
      if(size(state)/=nstate) error stop 'randtoolbox: MT state size mismatch'
      if(index<1 .or. index>nstate+2) error stop 'randtoolbox: MT state index invalid'
      this%state=modulo(state,4294967296_int64); this%idx=index
   end subroutine mt_put_state
end module randtoolbox_mt19937
