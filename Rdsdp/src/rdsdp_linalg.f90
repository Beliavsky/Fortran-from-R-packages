! LAPACK helpers used by the DSDP translation.
! DSDP copyright/license: see licenses/DSDP-LICENSE.
module rdsdp_linalg
   use rdsdp_kinds, only : dp
   implicit none
   private
   public :: spd_inverse_logdet, spd_logdet, solve_spd, frob_dot, norm2_dp, min_eigenvalue_sym

   interface
      subroutine dpotrf(uplo,n,a,lda,info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotrf
      subroutine dpotri(uplo,n,a,lda,info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotri
      subroutine dposv(uplo,n,nrhs,a,lda,b,ldb,info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: info
      end subroutine dposv
      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*), work(*)
         integer, intent(out) :: info
      end subroutine dsyev
   end interface

contains

   pure real(dp) function frob_dot(a,b) result(v)
      real(dp), intent(in) :: a(:,:), b(:,:)
      v = sum(a*b)
   end function frob_dot

   pure real(dp) function norm2_dp(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp,dot_product(x,x)))
   end function norm2_dp

   subroutine spd_inverse_logdet(a,ainv,logdet,ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      real(dp), intent(out) :: logdet
      logical, intent(out) :: ok
      real(dp), allocatable :: fac(:,:)
      integer :: n, i, j, info
      n = size(a,1)
      allocate(fac(n,n),ainv(n,n))
      fac = 0.5_dp*(a+transpose(a))
      call dpotrf('L',n,fac,n,info)
      if (info /= 0) then
         ok = .false.
         logdet = -huge(1.0_dp)
         ainv = 0.0_dp
         return
      end if
      logdet = 0.0_dp
      do i = 1, n
         logdet = logdet + 2.0_dp*log(fac(i,i))
      end do
      call dpotri('L',n,fac,n,info)
      if (info /= 0) then
         ok = .false.
         ainv = 0.0_dp
         return
      end if
      ainv = 0.0_dp
      do j = 1, n
         do i = j, n
            ainv(i,j) = fac(i,j)
            ainv(j,i) = fac(i,j)
         end do
      end do
      ok = .true.
   end subroutine spd_inverse_logdet

   subroutine spd_logdet(a,logdet,ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: logdet
      logical, intent(out) :: ok
      real(dp), allocatable :: fac(:,:)
      integer :: n,i,info
      n=size(a,1); allocate(fac(n,n)); fac=0.5_dp*(a+transpose(a))
      call dpotrf('L',n,fac,n,info)
      if (info/=0) then
         ok=.false.; logdet=-huge(1.0_dp); return
      end if
      logdet=0.0_dp
      do i=1,n; logdet=logdet+2.0_dp*log(fac(i,i)); end do
      ok=.true.
   end subroutine spd_logdet

   subroutine solve_spd(a,b,x,reg,ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: reg
      logical, intent(out) :: ok
      real(dp), allocatable :: aa(:,:), bb(:,:)
      real(dp) :: add
      integer :: n, i, info, attempt
      n = size(b)
      allocate(aa(n,n),bb(n,1))
      add = max(0.0_dp,reg)
      do attempt = 1, 9
         aa = 0.5_dp*(a+transpose(a))
         do i = 1, n
            aa(i,i) = aa(i,i) + add
         end do
         bb(:,1) = b
         call dposv('L',n,1,aa,n,bb,n,info)
         if (info == 0) then
            x = bb(:,1)
            ok = .true.
            return
         end if
         add = max(1.0e-14_dp,10.0_dp*max(add,1.0e-14_dp))
      end do
      x = 0.0_dp
      ok = .false.
   end subroutine solve_spd

   subroutine min_eigenvalue_sym(a,wmin,ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: wmin
      logical, intent(out) :: ok
      real(dp), allocatable :: aa(:,:), w(:), work(:)
      real(dp) :: qwork(1)
      integer :: n, lwork, info
      n = size(a,1)
      allocate(aa(n,n),w(n))
      aa = 0.5_dp*(a+transpose(a))
      lwork = -1
      call dsyev('N','L',n,aa,n,w,qwork,lwork,info)
      if (info /= 0) then
         ok = .false.; wmin = -huge(1.0_dp); return
      end if
      lwork = max(1,int(qwork(1)))
      allocate(work(lwork))
      aa = 0.5_dp*(a+transpose(a))
      call dsyev('N','L',n,aa,n,w,work,lwork,info)
      ok = info == 0
      if (ok) then
         wmin = w(1)
      else
         wmin = -huge(1.0_dp)
      end if
   end subroutine min_eigenvalue_sym

end module rdsdp_linalg
