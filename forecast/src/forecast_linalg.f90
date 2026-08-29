module forecast_linalg
   use forecast_kinds, only : dp
   implicit none
   private
   public :: least_squares, solve_linear, inverse_matrix, logdet_spd, covariance_ols, symmetric_eigen
   public :: polynomial_max_root_mod, cholesky_lower

   interface
      subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
         import dp
         character(len=1), intent(in) :: trans
         integer,intent(in) :: m,n,nrhs,lda,ldb,lwork
         real(dp),intent(inout) :: a(lda,*),b(ldb,*),work(*)
         integer,intent(out) :: info
      end subroutine
      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         import dp
         integer,intent(in) :: n,nrhs,lda,ldb
         real(dp),intent(inout) :: a(lda,*),b(ldb,*)
         integer,intent(out) :: ipiv(*),info
      end subroutine
      subroutine dpotrf(uplo,n,a,lda,info)
         import dp
         character(len=1),intent(in) :: uplo
         integer,intent(in) :: n,lda
         real(dp),intent(inout) :: a(lda,*)
         integer,intent(out) :: info
      end subroutine
      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import dp
         character(len=1),intent(in) :: jobz,uplo
         integer,intent(in) :: n,lda,lwork
         real(dp),intent(inout) :: a(lda,*),work(*)
         real(dp),intent(out) :: w(*)
         integer,intent(out) :: info
      end subroutine
      subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
         import dp
         character(len=1),intent(in) :: jobvl,jobvr
         integer,intent(in) :: n,lda,ldvl,ldvr,lwork
         real(dp),intent(inout) :: a(lda,*),work(*)
         real(dp),intent(out) :: wr(*),wi(*),vl(ldvl,*),vr(ldvr,*)
         integer,intent(out) :: info
      end subroutine
      subroutine dpotri(uplo,n,a,lda,info)
         import dp
         character(len=1),intent(in) :: uplo
         integer,intent(in) :: n,lda
         real(dp),intent(inout) :: a(lda,*)
         integer,intent(out) :: info
      end subroutine
   end interface
contains
   subroutine least_squares(x,y,beta,resid,rank,info)
      real(dp),intent(in) :: x(:,:),y(:)
      real(dp),allocatable,intent(out) :: beta(:),resid(:)
      integer,intent(out),optional :: rank,info
      real(dp),allocatable :: a(:,:),bb(:,:),work(:)
      real(dp) :: wk(1)
      integer :: m,n,ldb,lwork,istat
      m=size(x,1)
      n=size(x,2)
      ldb=max(m,n)
      allocate(a(m,n),bb(ldb,1))
      a=x
      bb=0.0_dp
      bb(1:m,1)=y
      call dgels('N',m,n,1,a,m,bb,ldb,wk,-1,istat)
      lwork=max(1,int(wk(1)))
      allocate(work(lwork))
      call dgels('N',m,n,1,a,m,bb,ldb,work,lwork,istat)
      allocate(beta(n),resid(m))
      beta=bb(1:n,1)
      resid=y-matmul(x,beta)
      if(present(rank)) rank=min(m,n)
      if(present(info)) info=istat
   end subroutine

   subroutine solve_linear(a,b,x,info)
      real(dp),intent(in) :: a(:,:),b(:)
      real(dp),allocatable,intent(out) :: x(:)
      integer,intent(out),optional :: info
      real(dp),allocatable :: aa(:,:),bb(:,:)
      integer,allocatable :: ipiv(:)
      integer :: n,istat
      n=size(a,1)
      allocate(aa(n,n),bb(n,1),ipiv(n))
      aa=a
      bb(:,1)=b
      call dgesv(n,1,aa,n,ipiv,bb,n,istat)
      allocate(x(n))
      x=bb(:,1)
      if(present(info)) info=istat
   end subroutine

   subroutine inverse_matrix(a,ainv,info)
      real(dp),intent(in) :: a(:,:)
      real(dp),allocatable,intent(out) :: ainv(:,:)
      integer,intent(out),optional :: info
      real(dp),allocatable :: aa(:,:),b(:,:),work(:)
      integer,allocatable :: ipiv(:)
      integer :: n,j,istat
      n=size(a,1)
      allocate(aa(n,n),b(n,n),ipiv(n))
      aa=a
      b=0.0_dp
      do j=1,n
      b(j,j)=1.0_dp
      end do
      call dgesv(n,n,aa,n,ipiv,b,n,istat)
      allocate(ainv(n,n))
      ainv=b
      if(present(info)) info=istat
   end subroutine

   subroutine cholesky_lower(a,l,info)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::l(:,:)
      integer,intent(out),optional::info
      integer::n,i,j,istat
      n=size(a,1)
      if(size(a,2)/=n)error stop 'cholesky_lower: matrix must be square'
      allocate(l(n,n))
      l=a
      call dpotrf('L',n,l,n,istat)
      if(istat==0)then
         do j=2,n
            do i=1,j-1
               l(i,j)=0.0_dp
            end do
         end do
      end if
      if(present(info))info=istat
   end subroutine cholesky_lower

   subroutine logdet_spd(a,logdet,info)
      real(dp),intent(in) :: a(:,:)
      real(dp),intent(out) :: logdet
      integer,intent(out),optional :: info
      real(dp),allocatable :: aa(:,:)
      integer :: n,j,istat
      n=size(a,1)
      allocate(aa(n,n))
      aa=a
      call dpotrf('L',n,aa,n,istat)
      logdet=0.0_dp
      if(istat==0) then
         do j=1,n
         logdet=logdet+2.0_dp*log(aa(j,j))
         end do
      else
         logdet=huge(1.0_dp)
      end if
      if(present(info)) info=istat
   end subroutine


   subroutine symmetric_eigen(a,values,vectors,info)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::values(:),vectors(:,:)
      integer,intent(out),optional::info
      real(dp),allocatable::aa(:,:),work(:)
      real(dp)::wk(1)
      integer::n,lwork,istat
      n=size(a,1)
      if(size(a,2)/=n)error stop 'symmetric_eigen: matrix must be square'
      allocate(aa(n,n),values(n))
      aa=a
      call dsyev('V','U',n,aa,n,values,wk,-1,istat)
      if(istat==0)then
         lwork=max(1,int(wk(1)))
         allocate(work(lwork))
         call dsyev('V','U',n,aa,n,values,work,lwork,istat)
      end if
      vectors=aa
      if(present(info))info=istat
   end subroutine symmetric_eigen


   function polynomial_max_root_mod(coeff,info) result(rmax)
      ! coeff contains c0,...,cn for c0+c1*z+...+cn*z**n.
      real(dp),intent(in)::coeff(:)
      integer,intent(out),optional::info
      real(dp)::rmax,wk(1),dummy_l(1,1),dummy_r(1,1)
      real(dp),allocatable::a(:,:),wr(:),wi(:),work(:)
      integer::n,j,lwork,istat
      n=size(coeff)-1
      rmax=0.0_dp
      if(n<1 .or. abs(coeff(size(coeff)))<=tiny(1.0_dp))then
         if(present(info))info=1
         rmax=huge(1.0_dp)
         return
      end if
      allocate(a(n,n),wr(n),wi(n))
      a=0.0_dp
      do j=1,n
         a(1,j)=-coeff(n-j+1)/coeff(n+1)
      end do
      do j=2,n
         a(j,j-1)=1.0_dp
      end do
      call dgeev('N','N',n,a,n,wr,wi,dummy_l,1,dummy_r,1,wk,-1,istat)
      if(istat==0)then
         lwork=max(1,int(wk(1)))
         allocate(work(lwork))
         ! dgeev overwrites A, so rebuild after workspace query.
         a=0.0_dp
         do j=1,n
            a(1,j)=-coeff(n-j+1)/coeff(n+1)
         end do
         do j=2,n
            a(j,j-1)=1.0_dp
         end do
         call dgeev('N','N',n,a,n,wr,wi,dummy_l,1,dummy_r,1,work,lwork,istat)
      end if
      if(istat==0)then
         rmax=maxval(sqrt(wr*wr+wi*wi))
      else
         rmax=huge(1.0_dp)
      end if
      if(present(info))info=istat
   end function polynomial_max_root_mod

   subroutine covariance_ols(x,sigma2,cov,info)
      real(dp),intent(in) :: x(:,:),sigma2
      real(dp),allocatable,intent(out) :: cov(:,:)
      integer,intent(out),optional :: info
      real(dp),allocatable :: xtx(:,:),inv(:,:)
      integer :: istat
      xtx=matmul(transpose(x),x)
      call inverse_matrix(xtx,inv,istat)
      cov=sigma2*inv
      if(present(info)) info=istat
   end subroutine
end module forecast_linalg
