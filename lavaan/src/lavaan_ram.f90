module lavaan_ram
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general
   implicit none
   private
   integer, parameter, public :: ram_a=1, ram_s=2, ram_m=3
   type, public :: ram_model
      real(dp), allocatable :: a(:, :), s(:, :), m(:)
      integer, allocatable :: observed(:)
   end type ram_model
   type, public :: ram_free_map
      integer, allocatable :: matrix_id(:), row(:), col(:)
   end type ram_free_map
   public :: ram_sigma, ram_mu, ram_set_free, ram_get_free, ram_validate
contains
   subroutine ram_validate(model, ok)
      type(ram_model), intent(in)::model
      logical,intent(out)::ok
      integer::n
      ok=.false.
      if(.not.allocated(model%a) .or. .not.allocated(model%s)) return
      n=size(model%a,1)
      if(size(model%a,2)/=n .or. any(shape(model%s)/=[n,n])) return
      if(.not.allocated(model%observed)) return
      if(any(model%observed<1) .or. any(model%observed>n)) return
      ok=.true.
   end subroutine ram_validate

   subroutine ram_sigma(model, sigma, info)
      type(ram_model), intent(in)::model
      real(dp),allocatable,intent(out)::sigma(:,:)
      integer,intent(out)::info
      real(dp),allocatable::ia(:,:),inv(:,:),allcov(:,:)
      integer::n,p,i,j
      logical::ok
      call ram_validate(model,ok)
      if(.not.ok) then
      info=-1
      allocate(sigma(0,0))
      return
      end if
      n=size(model%a,1)
      p=size(model%observed)
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1
      end do
      call inverse_general(ia,inv,info)
      if(info/=0) then
      allocate(sigma(p,p))
      sigma=huge(1.0_dp)
      return
      end if
      allcov=matmul(inv,matmul(model%s,transpose(inv)))
      allocate(sigma(p,p))
      do j=1,p
      do i=1,p
      sigma(i,j)=allcov(model%observed(i),model%observed(j))
      end do
      end do
   end subroutine ram_sigma

   subroutine ram_mu(model, mu, info)
      type(ram_model), intent(in)::model
      real(dp),allocatable,intent(out)::mu(:)
      integer,intent(out)::info
      real(dp),allocatable::ia(:,:),inv(:,:),allmu(:)
      integer::n,p,i
      logical::ok
      call ram_validate(model,ok)
      if(.not.ok) then
      info=-1
      allocate(mu(0))
      return
      end if
      n=size(model%a,1)
      p=size(model%observed)
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1
      end do
      call inverse_general(ia,inv,info)
      allocate(mu(p))
      if(info/=0) then
      mu=huge(1.0_dp)
      return
      end if
      if(allocated(model%m)) then
      allmu=matmul(inv,model%m)
      else
      allocate(allmu(n))
      allmu=0
      end if
      do i=1,p
      mu(i)=allmu(model%observed(i))
      end do
   end subroutine ram_mu

   function ram_get_free(model,map) result(x)
      type(ram_model),intent(in)::model
      type(ram_free_map),intent(in)::map
      real(dp),allocatable::x(:)
      integer::k
      allocate(x(size(map%matrix_id)))
      do k=1,size(x)
         select case(map%matrix_id(k))
         case(ram_a); x(k)=model%a(map%row(k),map%col(k))
         case(ram_s); x(k)=model%s(map%row(k),map%col(k))
         case(ram_m); x(k)=model%m(map%row(k))
         end select
      end do
   end function ram_get_free

   subroutine ram_set_free(model,map,x)
      type(ram_model),intent(inout)::model
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::x(:)
      integer::k,r,c
      do k=1,size(x)
         r=map%row(k)
         c=map%col(k)
         select case(map%matrix_id(k))
         case(ram_a); model%a(r,c)=x(k)
         case(ram_s)
         model%s(r,c)=x(k)
         model%s(c,r)=x(k)
         case(ram_m); model%m(r)=x(k)
         end select
      end do
   end subroutine ram_set_free
end module lavaan_ram
