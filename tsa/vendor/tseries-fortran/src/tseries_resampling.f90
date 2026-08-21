! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_resampling
   use tseries_kinds, only : dp
   use tseries_random, only : random_integer, random_exponential, random_permutation, random_normal
   use tseries_fft, only : discrete_fourier_transform, inverse_discrete_fourier_transform
   use tseries_stats, only : rank_values
   implicit none
   private

   public :: quadratic_map
   public :: stationary_bootstrap
   public :: block_bootstrap
   public :: permutation_surrogate
   public :: fft_surrogate
   public :: amplitude_surrogate

contains

   function quadratic_map(xi,a,n) result(x)
      real(dp), intent(in) :: xi,a
      integer, intent(in) :: n
      real(dp), allocatable :: x(:)
      integer :: i
      allocate(x(max(0,n)))
      if(n<=0) return
      x(1)=xi
      do i=2,n
         x(i)=a*(1.0_dp-x(i-1))*x(i-1)
      end do
   end function quadratic_map

   subroutine stationary_bootstrap(x,mean_block_length,sample)
      real(dp), intent(in) :: x(:),mean_block_length
      real(dp), intent(out) :: sample(:)
      real(dp) :: p
      integer :: n,i,start_index,block_length,j,index
      n=size(x)
      if(size(sample)/=n .or. n==0) return
      p=1.0_dp/max(mean_block_length,1.0_dp)
      i=1
      do while(i<=n)
         start_index=random_integer(1,n)
         block_length=1+int(random_exponential()/max(-log(max(1.0e-12_dp,1.0_dp-p)),1.0e-12_dp))
         do j=0,block_length-1
            if(i>n) exit
            index=1+modulo(start_index-1+j,n)
            sample(i)=x(index)
            i=i+1
         end do
      end do
   end subroutine stationary_bootstrap

   subroutine block_bootstrap(x,block_length,sample)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: block_length
      real(dp), intent(out) :: sample(:)
      integer :: n,i,start_index,j,l
      n=size(x); l=max(1,min(block_length,n))
      if(size(sample)/=n .or. n==0) return
      i=1
      do while(i<=n)
         start_index=random_integer(1,n-l+1)
         do j=0,l-1
            if(i>n) exit
            sample(i)=x(start_index+j)
            i=i+1
         end do
      end do
   end subroutine block_bootstrap

   subroutine permutation_surrogate(x,y)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
      integer, allocatable :: p(:)
      integer :: n
      n=size(x)
      if(size(y)/=n) return
      allocate(p(n)); call random_permutation(n,p)
      y=x(p)
   end subroutine permutation_surrogate

   subroutine fft_surrogate(x,y)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
      complex(dp), allocatable :: z(:)
      complex(dp) :: phase
      real(dp) :: angle
      integer :: n,k,pair
      real(dp), parameter :: two_pi=2.0_dp*acos(-1.0_dp)
      n=size(x)
      if(size(y)/=n) return
      allocate(z(n)); call discrete_fourier_transform(x,z)
      do k=2,(n+1)/2
         call random_number(angle); angle=two_pi*angle
         phase=cmplx(cos(angle),sin(angle),dp)
         z(k)=z(k)*phase
         pair=n-k+2
         if(pair>=1 .and. pair<=n .and. pair/=k) z(pair)=conjg(z(k))
      end do
      if(mod(n,2)==0) z(n/2+1)=cmplx(real(z(n/2+1),dp),0.0_dp,dp)
      call inverse_discrete_fourier_transform(z,y)
   end subroutine fft_surrogate

   subroutine amplitude_surrogate(x,y)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
      real(dp), allocatable :: sorted_x(:),g(:),sorted_g(:),gaussianized(:),phase_random(:)
      integer, allocatable :: ranks_x(:),ranks_y(:)
      integer :: n,i
      n=size(x)
      if(size(y)/=n) return
      allocate(sorted_x(n),g(n),sorted_g(n),gaussianized(n),phase_random(n),ranks_x(n),ranks_y(n))
      sorted_x=x; call sort_real(sorted_x)
      do i=1,n
         g(i)=random_normal()
      end do
      sorted_g=g; call sort_real(sorted_g)
      call rank_values(x,ranks_x)
      do i=1,n
         gaussianized(i)=sorted_g(ranks_x(i))
      end do
      call fft_surrogate(gaussianized,phase_random)
      call rank_values(phase_random,ranks_y)
      do i=1,n
         y(i)=sorted_x(ranks_y(i))
      end do
   end subroutine amplitude_surrogate

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i,j
      do i=2,size(x)
         key=x(i); j=i-1
         do while(j>=1)
            if(x(j)<=key) exit
            x(j+1)=x(j); j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real

end module tseries_resampling
