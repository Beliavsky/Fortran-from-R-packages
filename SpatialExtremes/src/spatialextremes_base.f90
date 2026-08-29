module spatialextremes_base
   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_compat, only: rexp, rchisq
   use r_linalg, only: shared_cholesky_factor => cholesky_factor
   use r_linalg, only: shared_solve_spd => solve_spd
   use r_linalg, only: shared_spd_inverse_logdet => spd_inverse_logdet
   implicit none
   private
   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: neg_huge = -1.0e15_dp
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   public :: pair_count, pair_indices, euclidean_distances, distance_vectors, distance_to_origin
   public :: chol_upper, solve_spd, inverse_spd, logdet_spd, is_finite, nan_dp, eye, exp_rand, chisq_rand

contains
   pure integer function pair_count(n) result(np)
      integer, intent(in) :: n
      np = n*(n-1)/2
   end function pair_count

   pure subroutine pair_indices(k,n,i,j)
      integer, intent(in) :: k,n
      integer, intent(out) :: i,j
      integer :: c, ii, jj
      c=0
      do ii=1,n-1
         do jj=ii+1,n
            c=c+1
            if (c==k) then
               i=ii
               j=jj
               return
            end if
         end do
      end do
      i=0
      j=0
   end subroutine pair_indices

   function euclidean_distances(coord) result(d)
      real(dp), intent(in) :: coord(:,:)
      real(dp), allocatable :: d(:)
      integer :: n,np,k,i,j
      n=size(coord,1)
      np=pair_count(n)
      allocate(d(np))
      do k=1,np
         call pair_indices(k,n,i,j)
         d(k)=sqrt(sum((coord(j,:)-coord(i,:))**2))
      end do
   end function euclidean_distances

   function distance_vectors(coord) result(dv)
      real(dp), intent(in) :: coord(:,:)
      real(dp), allocatable :: dv(:,:)
      integer :: n,np,k,i,j
      n=size(coord,1)
      np=pair_count(n)
      allocate(dv(np,size(coord,2)))
      do k=1,np
         call pair_indices(k,n,i,j)
         dv(k,:)=coord(j,:)-coord(i,:)
      end do
   end function distance_vectors

   function distance_to_origin(coord) result(d)
      real(dp), intent(in) :: coord(:,:)
      real(dp), allocatable :: d(:)
      integer :: i
      allocate(d(size(coord,1)))
      do i=1,size(coord,1)
         d(i)=sqrt(sum(coord(i,:)**2))
      end do
   end function distance_to_origin

   logical pure function is_finite(x) result(ok)
      real(dp), intent(in) :: x
      ok=ieee_is_finite(x)
   end function is_finite

   pure real(dp) function nan_dp() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nan_dp

   real(dp) function exp_rand() result(x)
      real(dp), allocatable :: a(:)
      a=rexp(1,1.0_dp)
      x=a(1)
   end function exp_rand

   real(dp) function chisq_rand(df) result(x)
      real(dp),intent(in)::df
      real(dp), allocatable :: a(:)
      a=rchisq(1,df)
      x=a(1)
   end function chisq_rand

   function eye(n) result(a)
      integer,intent(in)::n
      real(dp)::a(n,n)
      integer::i
      a=0.0_dp
      do i=1,n
      a(i,i)=1.0_dp
      end do
   end function eye

   subroutine chol_upper(a,r,info)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::r(size(a,1),size(a,2))
      integer,intent(out)::info
      real(dp), allocatable :: factor(:,:)
      call shared_cholesky_factor(a,factor,info,upper=.true.)
      if(info==0) r=factor
   end subroutine chol_upper

   subroutine solve_spd(a,b,x,info)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp),intent(out)::x(size(b,1),size(b,2))
      integer,intent(out)::info
      call shared_solve_spd(a,b,x,info,upper=.true.)
   end subroutine solve_spd

   subroutine inverse_spd(a,ainv,info)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::ainv(size(a,1),size(a,2))
      integer,intent(out)::info
      real(dp), allocatable :: inverse(:,:)
      real(dp) :: ignored_logdet
      call shared_spd_inverse_logdet(a,inverse,ignored_logdet,info)
      if(info==0) ainv=inverse
   end subroutine inverse_spd

   real(dp) function logdet_spd(a,info) result(v)
      real(dp),intent(in)::a(:,:)
      integer,intent(out)::info
      real(dp), allocatable :: inverse(:,:)
      call shared_spd_inverse_logdet(a,inverse,v,info)
      if(info/=0) v=nan_dp()
   end function logdet_spd
end module spatialextremes_base
