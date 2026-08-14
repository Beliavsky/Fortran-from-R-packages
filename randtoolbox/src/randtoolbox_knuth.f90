! Knuth's subtractive/additive lagged-Fibonacci generator from TAOCP 3.6.
! The upstream C code states: public domain and freely copyable.
module randtoolbox_knuth
   use, intrinsic :: iso_fortran_env, only : int64, real64
   implicit none
   private
   integer, parameter :: kk=100, ll=37

   type, public :: knuth_rng
      real(real64) :: state(kk)=0.0_real64
      logical :: initialized=.false.
   contains
      procedure, public :: seed => knuth_seed
      procedure, public :: fill => knuth_fill
      procedure, public :: fill_matrix => knuth_fill_matrix
   end type knuth_rng

contains
   pure real(real64) function mod_sum(x,y) result(z)
      real(real64), intent(in) :: x,y
      z = (x+y) - real(int(x+y),real64)
   end function mod_sum

   subroutine knuth_array(this,a)
      class(knuth_rng), intent(inout) :: this
      real(real64), intent(out) :: a(:)
      integer :: i,j,n
      n=size(a)
      if(n < kk) error stop 'randtoolbox: Knuth array length must be at least 100'
      a(1:kk)=this%state
      do j=kk+1,n
         a(j)=mod_sum(a(j-kk),a(j-ll))
      end do
      j=n+1
      do i=1,ll
         this%state(i)=mod_sum(a(j-kk),a(j-ll)); j=j+1
      end do
      do i=ll+1,kk
         this%state(i)=mod_sum(a(j-kk),this%state(i-ll)); j=j+1
      end do
   end subroutine knuth_array

   subroutine knuth_seed(this,seed)
      class(knuth_rng), intent(inout) :: this
      integer(int64), intent(in) :: seed
      integer :: t,s,j
      real(real64) :: u(2*kk-1), ulp, ss
      integer(int64) :: smask
      ulp=(1.0_real64/real(2_int64**30,real64))/real(2_int64**22,real64)
      smask=iand(seed,int(z'3fffffff',int64))
      ss=2.0_real64*ulp*real(smask+2_int64,real64)
      do j=1,kk
         u(j)=ss
         ss=ss+ss
         if(ss>=1.0_real64) ss=ss-(1.0_real64-2.0_real64*ulp)
      end do
      u(2)=u(2)+ulp
      s=int(smask)
      t=69
      do while(t/=0)
         do j=kk,2,-1
            u(2*j-1)=u(j)
            u(2*j-2)=0.0_real64
         end do
         do j=2*kk-1,kk+1,-1
            u(j-(kk-ll))=mod_sum(u(j-(kk-ll)),u(j))
            u(j-kk)=mod_sum(u(j-kk),u(j))
         end do
         if(iand(s,1)/=0) then
            do j=kk+1,2,-1
               u(j)=u(j-1)
            end do
            u(1)=u(kk+1)
            u(ll+1)=mod_sum(u(ll+1),u(kk+1))
         end if
         if(s/=0) then
            s=shiftr(s,1)
         else
            t=t-1
         end if
      end do
      do j=1,ll
         this%state(j+kk-ll)=u(j)
      end do
      do j=ll+1,kk
         this%state(j-ll)=u(j)
      end do
      do j=1,10
         call knuth_array(this,u)
      end do
      this%initialized=.true.
   end subroutine knuth_seed

   subroutine knuth_fill(this,x)
      class(knuth_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:)
      real(real64), allocatable :: tmp(:)
      integer :: n
      if(.not.this%initialized) call this%seed(314159_int64)
      n=size(x)
      if(n==0) return
      if(n<=kk) then
         allocate(tmp(kk+1)); call knuth_array(this,tmp); x=tmp(1:n)
      else
         call knuth_array(this,x)
      end if
   end subroutine knuth_fill

   subroutine knuth_fill_matrix(this,x)
      class(knuth_rng), intent(inout) :: this
      real(real64), intent(out) :: x(:,:)
      real(real64), allocatable :: v(:)
      integer :: i,j,k
      allocate(v(size(x))); call this%fill(v)
      k=0
      do j=1,size(x,2)
         do i=1,size(x,1)
            k=k+1; x(i,j)=v(k)
         end do
      end do
   end subroutine knuth_fill_matrix
end module randtoolbox_knuth
