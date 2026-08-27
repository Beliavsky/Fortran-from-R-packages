! SPDX-License-Identifier: BSD-3-Clause
! External-name adapters used by the fixed-form ARPACK sources.
subroutine daxpy(n,a,x,ix,y,iy)
  use la_blas_d, only: backend => la_daxpy
  integer :: n,ix,iy
  double precision :: a,x(*),y(*)
  call backend(n,a,x,ix,y,iy)
end subroutine daxpy

subroutine dcopy(n,x,ix,y,iy)
  use la_blas_d, only: backend => la_dcopy
  integer :: n,ix,iy
  double precision :: x(*),y(*)
  call backend(n,x,ix,y,iy)
end subroutine dcopy

double precision function ddot(n,x,ix,y,iy)
  use la_blas_d, only: backend => la_ddot
  integer :: n,ix,iy
  double precision :: x(*),y(*)
  ddot=backend(n,x,ix,y,iy)
end function ddot

subroutine dgemv(t,m,n,a,x,ldx,y,iy,b,z,iz)
  use la_blas_d, only: backend => la_dgemv
  character(len=1) :: t
  integer :: m,n,ldx,iy,iz
  double precision :: a,b,x(ldx,*),y(*),z(*)
  call backend(t,m,n,a,x,ldx,y,iy,b,z,iz)
end subroutine dgemv

subroutine dger(m,n,a,x,ix,y,iy,z,ldz)
  use la_blas_d, only: backend => la_dger
  integer :: m,n,ix,iy,ldz
  double precision :: a,x(*),y(*),z(ldz,*)
  call backend(m,n,a,x,ix,y,iy,z,ldz)
end subroutine dger

double precision function dnrm2(n,x,ix)
  use la_blas_d, only: backend => la_dnrm2
  integer :: n,ix
  double precision :: x(*)
  dnrm2=backend(n,x,ix)
end function dnrm2

subroutine drot(n,x,ix,y,iy,c,s)
  use la_blas_d, only: backend => la_drot
  integer :: n,ix,iy
  double precision :: x(*),y(*),c,s
  call backend(n,x,ix,y,iy,c,s)
end subroutine drot

subroutine dscal(n,a,x,ix)
  use la_blas_d, only: backend => la_dscal
  integer :: n,ix
  double precision :: a,x(*)
  call backend(n,a,x,ix)
end subroutine dscal

subroutine dswap(n,x,ix,y,iy)
  use la_blas_d, only: backend => la_dswap
  integer :: n,ix,iy
  double precision :: x(*),y(*)
  call backend(n,x,ix,y,iy)
end subroutine dswap

subroutine dtrmm(si,u,t,d,m,n,a,x,ldx,y,ldy)
  use la_blas_d, only: backend => la_dtrmm
  character(len=1) :: si,u,t,d
  integer :: m,n,ldx,ldy
  double precision :: a,x(ldx,*),y(ldy,*)
  call backend(si,u,t,d,m,n,a,x,ldx,y,ldy)
end subroutine dtrmm

subroutine dlabad(small,large)
  use la_lapack_d, only: backend => la_dlabad
  double precision :: small,large
  call backend(small,large)
end subroutine dlabad

subroutine dlacpy(u,m,n,a,lda,b,ldb)
  use la_lapack_d, only: backend => la_dlacpy
  character(len=1) :: u
  integer :: m,n,lda,ldb
  double precision :: a(lda,*),b(ldb,*)
  call backend(u,m,n,a,lda,b,ldb)
end subroutine dlacpy

subroutine dlae2(a,b,c,r1,r2)
  use la_lapack_d, only: backend => la_dlae2
  double precision :: a,b,c,r1,r2
  call backend(a,b,c,r1,r2)
end subroutine dlae2

subroutine dlaev2(a,b,c,r1,r2,cs,sn)
  use la_lapack_d, only: backend => la_dlaev2
  double precision :: a,b,c,r1,r2,cs,sn
  call backend(a,b,c,r1,r2,cs,sn)
end subroutine dlaev2

subroutine dgeqr2(m,n,a,lda,tau,work,info)
  use la_lapack_d, only: backend => la_dgeqr2
  integer :: m,n,lda,info
  double precision :: a(lda,*),tau(*),work(*)
  call backend(m,n,a,lda,tau,work,info)
end subroutine dgeqr2

subroutine dlahqr(wt,wz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,info)
  use la_lapack_d, only: backend => la_dlahqr
  logical :: wt,wz
  integer :: n,ilo,ihi,ldh,iloz,ihiz,ldz,info
  double precision :: h(ldh,*),wr(*),wi(*),z(ldz,*)
  call backend(wt,wz,n,ilo,ihi,h,ldh,wr,wi,iloz,ihiz,z,ldz,info)
end subroutine dlahqr

double precision function dlamch(c)
  use la_lapack_d, only: backend => la_dlamch
  character(len=*) :: c
  dlamch=backend(c)
end function dlamch

double precision function dlanhs(nm,n,a,lda,work)
  use la_lapack_d, only: backend => la_dlanhs
  character(len=1) :: nm
  integer :: n,lda
  double precision :: a(lda,*),work(*)
  dlanhs=backend(nm,n,a,lda,work)
end function dlanhs

double precision function dlanst(nm,n,d,e)
  use la_lapack_d, only: backend => la_dlanst
  character(len=1) :: nm
  integer :: n
  double precision :: d(*),e(*)
  dlanst=backend(nm,n,d,e)
end function dlanst

subroutine dlanv2(a,b,c,d,r1r,r1i,r2r,r2i,cs,sn)
  use la_lapack_d, only: backend => la_dlanv2
  double precision :: a,b,c,d,r1r,r1i,r2r,r2i,cs,sn
  call backend(a,b,c,d,r1r,r1i,r2r,r2i,cs,sn)
end subroutine dlanv2

double precision function dlapy2(x,y)
  use la_lapack_d, only: backend => la_dlapy2
  double precision :: x,y
  dlapy2=backend(x,y)
end function dlapy2

subroutine dlarf(si,m,n,v,iv,tau,c,ldc,work)
  use la_lapack_d, only: backend => la_dlarf
  character(len=1) :: si
  integer :: m,n,iv,ldc
  double precision :: v(*),tau,c(ldc,*),work(*)
  call backend(si,m,n,v,iv,tau,c,ldc,work)
end subroutine dlarf

subroutine dlarfg(n,a,x,ix,tau)
  use la_lapack_d, only: backend => la_dlarfg
  integer :: n,ix
  double precision :: a,x(*),tau
  call backend(n,a,x,ix,tau)
end subroutine dlarfg

subroutine dlarnv(id,seed,n,x)
  use la_lapack_d, only: backend => la_dlarnv
  integer :: id,seed(*),n
  double precision :: x(*)
  call backend(id,seed,n,x)
end subroutine dlarnv

subroutine dlartg(f,g,c,s,r)
  use la_lapack_d, only: backend => la_dlartg
  double precision :: f,g,c,s,r
  call backend(f,g,c,s,r)
end subroutine dlartg

subroutine dlascl(ty,kl,ku,cf,ct,m,n,a,lda,info)
  use la_lapack_d, only: backend => la_dlascl
  character(len=1) :: ty
  integer :: kl,ku,m,n,lda,info
  double precision :: cf,ct,a(lda,*)
  call backend(ty,kl,ku,cf,ct,m,n,a,lda,info)
end subroutine dlascl

subroutine dlaset(u,m,n,a,b,x,ldx)
  use la_lapack_d, only: backend => la_dlaset
  character(len=1) :: u
  integer :: m,n,ldx
  double precision :: a,b,x(ldx,*)
  call backend(u,m,n,a,b,x,ldx)
end subroutine dlaset

subroutine dlasr(si,p,d,m,n,c,s,a,lda)
  use la_lapack_d, only: backend => la_dlasr
  character(len=1) :: si,p,d
  integer :: m,n,lda
  double precision :: c(*),s(*),a(lda,*)
  call backend(si,p,d,m,n,c,s,a,lda)
end subroutine dlasr

subroutine dlasrt(id,n,d,info)
  use la_lapack_d, only: backend => la_dlasrt
  character(len=1) :: id
  integer :: n,info
  double precision :: d(*)
  call backend(id,n,d,info)
end subroutine dlasrt

subroutine dorm2r(si,tr,m,n,k,a,lda,tau,c,ldc,work,info)
  use la_lapack_d, only: backend => la_dorm2r
  character(len=1) :: si,tr
  integer :: m,n,k,lda,ldc,info
  double precision :: a(lda,*),tau(*),c(ldc,*),work(*)
  call backend(si,tr,m,n,k,a,lda,tau,c,ldc,work,info)
end subroutine dorm2r

subroutine dsteqr(cz,n,d,e,z,ldz,work,info)
  use la_lapack_d, only: backend => la_dsteqr
  character(len=1) :: cz
  integer :: n,ldz,info
  double precision :: d(*),e(*),z(ldz,*),work(*)
  call backend(cz,n,d,e,z,ldz,work,info)
end subroutine dsteqr

subroutine dtrevc(si,hm,sel,n,t,ldt,vl,ldvl,vr,ldvr,mm,m,work,info)
  use la_lapack_d, only: backend => la_dtrevc
  character(len=1) :: si,hm
  logical :: sel(*)
  integer :: n,ldt,ldvl,ldvr,mm,m,info
  double precision :: t(ldt,*),vl(ldvl,*),vr(ldvr,*),work(*)
  call backend(si,hm,sel,n,t,ldt,vl,ldvl,vr,ldvr,mm,m,work,info)
end subroutine dtrevc

subroutine dtrsen(j,cq,sel,n,t,ldt,q,ldq,wr,wi,m,s,sep,w,lw,iw,liw,info)
  use la_lapack_d, only: backend => la_dtrsen
  character(len=1) :: j,cq
  logical :: sel(*)
  integer :: n,ldt,ldq,m,lw,iw(*),liw,info
  double precision :: t(ldt,*),q(ldq,*),wr(*),wi(*),s,sep,w(*)
  call backend(j,cq,sel,n,t,ldt,q,ldq,wr,wi,m,s,sep,w,lw,iw,liw,info)
end subroutine dtrsen

logical function lsame(a,b)
  use la_blas_aux, only: backend => la_lsame
  character(len=1) :: a,b
  lsame=backend(a,b)
end function lsame

subroutine xerbla(name,info)
  use la_blas_aux, only: backend => la_xerbla
  character(len=*) :: name
  integer :: info
  call backend(name,info)
end subroutine xerbla
