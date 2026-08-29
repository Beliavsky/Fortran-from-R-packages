module lavaan_simulation
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_sigma, ram_mu
   use lavaan_linalg, only : chol_lower
   implicit none
   private
   public :: simulate_ram, random_seed_lavaan
contains
   subroutine random_seed_lavaan(seed)
      integer,intent(in)::seed
      integer::n,i
      integer,allocatable::put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i=1,n
      put(i)=seed+104729*i
      end do
      call random_seed(put=put)
   end subroutine random_seed_lavaan
   function normal_random() result(z)
      real(dp)::z,u1,u2
      call random_number(u1)
      call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2*log(u1))*cos(2*acos(-1.0_dp)*u2)
   end function normal_random
   subroutine simulate_ram(model,n,x,info)
      type(ram_model),intent(in)::model
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::x(:,:)
      integer,intent(out)::info
      real(dp),allocatable::sigma(:,:),mu(:),l(:,:),z(:)
      integer::i,j,p
      call ram_sigma(model,sigma,info)
      if(info/=0) then
      allocate(x(0,0))
      return
      end if
      call ram_mu(model,mu,info)
      call chol_lower(sigma,l,info)
      p=size(mu)
      allocate(x(n,p),z(p))
      if(info/=0) return
      do i=1,n
      do j=1,p
      z(j)=normal_random()
      end do
      x(i,:)=mu+matmul(l,z)
      end do
   end subroutine simulate_ram
end module lavaan_simulation
