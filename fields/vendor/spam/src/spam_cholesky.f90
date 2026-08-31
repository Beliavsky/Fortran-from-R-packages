module spam_cholesky
use spam_kinds, only: dp
use spam_types, only: csr_matrix, spam_chol
implicit none
private
public :: spam_chol_factor, spam_chol_update, spam_solve, spam_forwardsolve, spam_backsolve
public :: spam_logdet, spam_determinant, chol_factor_logdet, chol_to_dense, spam_chol2inv
interface
 subroutine cholstepwise(m,nnzd,d,jd,id,doperm,invp,perm,nsub,nsubmax,lindx,xlindx,nsuper,nnzlmax,lnz,xlnz,snode,xsuper,cachsz,ierr)
  integer m,nnzd,doperm,nsub,nsubmax,nsuper,nnzlmax,cachsz,ierr
  integer jd(nnzd),id(m+1),invp(m),perm(m),lindx(*),xlindx(*),xlnz(m+1),snode(m),xsuper(m+1)
  double precision d(nnzd),lnz(*)
 end subroutine
 subroutine updatefactor(m,nnzd,d,jd,id,invp,perm,lindx,xlindx,nsuper,lnz,xlnz,snode,xsuper,cachesize,ierr)
  integer m,nnzd,nsuper,ierr,jd(nnzd),cachesize,id(m+1),lindx(*),xlindx(*),invp(m),perm(m),xlnz(m+1),snode(m),xsuper(m+1)
  double precision d(nnzd),lnz(*)
 end subroutine
 subroutine backsolves(m,nsuper,nrhs,lindx,xlindx,lnz,xlnz,invp,perm,xsuper,newrhs,sol,b)
  integer m,nsuper,nrhs,lindx(*),xlindx(m+1),invp(m),perm(m),xlnz(m+1),xsuper(m+1)
  double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)
 end subroutine
 subroutine pivotforwardsolve(m,nsuper,nrhs,lindx,xlindx,lnz,xlnz,invp,perm,xsuper,newrhs,sol,b)
  integer m,nsuper,nrhs,lindx(*),xlindx(m+1),invp(m),perm(m),xlnz(m+1),xsuper(m+1)
  double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)
 end subroutine
 subroutine pivotbacksolve(m,nsuper,nrhs,lindx,xlindx,lnz,xlnz,invp,perm,xsuper,newrhs,sol,b)
  integer m,nsuper,nrhs,lindx(*),xlindx(m+1),invp(m),perm(m),xlnz(m+1),xsuper(m+1)
  double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)
 end subroutine
end interface
interface spam_solve
 module procedure spam_solve_vec, spam_solve_mat
end interface
contains

function spam_chol_factor(a,pivot,cache_kb) result(f)
type(csr_matrix),intent(in)::a
character(len=*),intent(in),optional::pivot
integer,intent(in),optional::cache_kb
type(spam_chol)::f
integer::n,nnz,doperm,nsub,nsubmax,nnzlmax,nsuper,ierr,cache,tries
integer,allocatable::invp(:),perm(:),lindx(:),xlindx(:),xlnz(:),snode(:),xsuper(:)
real(dp),allocatable::lnz(:)
character(len=16)::piv
if(a%nrow/=a%ncol) error stop 'spam_chol_factor: matrix must be square'
if(.not.a%valid()) error stop 'spam_chol_factor: invalid CSR matrix'
n=a%nrow;nnz=a%nnz();cache=512;if(present(cache_kb))cache=cache_kb
piv='mmd';if(present(pivot))piv=lower(trim(pivot))
select case(trim(piv));case('none','false','identity');doperm=0;case('rcm');doperm=2;case default;doperm=1;end select
nsubmax=max(n+1,4*max(nnz,n));nnzlmax=max(n+1,6*max(nnz,n));ierr=0
allocate(invp(n),perm(n),xlindx(n+1),xlnz(n+1),snode(n),xsuper(n+1))
if (doperm == 0) then
   perm = [(tries, tries=1,n)]
   invp = perm
end if
do tries=1,8
 allocate(lindx(nsubmax),lnz(nnzlmax));lindx=0;lnz=0.0_dp
 call cholstepwise(n,nnz,a%entries,a%colindices,a%rowpointers,doperm,invp,perm, &
      nsub,nsubmax,lindx,xlindx,nsuper,nnzlmax,lnz,xlnz,snode,xsuper,cache,ierr)
 if(ierr==0)exit
 deallocate(lindx,lnz)
 if(ierr==4)then;nnzlmax=2*nnzlmax;else if(ierr==5)then;nsubmax=2*nsubmax;else;exit;end if
end do
f%n=n;f%nsuper=max(0,nsuper);f%nnz_a=nnz;f%cache_kb=cache;f%info=ierr
if(ierr/=0)return
allocate(f%entries(xlnz(n+1)-1),f%rowpointers(n+1),f%pivot(n),f%invpivot(n),f%snmember(n))
allocate(f%colindices(nsub),f%colpointers(n+1),f%supernodes(n+1))
f%entries=lnz(:size(f%entries));f%rowpointers=xlnz;f%pivot=perm;f%invpivot=invp;f%snmember=snode
f%colpointers=0; f%supernodes=0
f%colindices=lindx(:nsub);f%colpointers(:nsuper+1)=xlindx(:nsuper+1);f%supernodes(:nsuper+1)=xsuper(:nsuper+1)
end function spam_chol_factor

subroutine spam_chol_update(f,a)
type(spam_chol),intent(inout)::f
type(csr_matrix),intent(in)::a
integer::ierr
if(f%info/=0) error stop 'spam_chol_update: invalid factor'
if(a%nrow/=f%n .or. a%ncol/=f%n) error stop 'spam_chol_update: dimension mismatch'
call updatefactor(f%n,a%nnz(),a%entries,a%colindices,a%rowpointers,f%invpivot,f%pivot,f%colindices,f%colpointers, &
 f%nsuper,f%entries,f%rowpointers,f%snmember,f%supernodes,f%cache_kb,ierr)
f%info=ierr
end subroutine spam_chol_update

function spam_solve_vec(f,b) result(x)
type(spam_chol),intent(in)::f;real(dp),intent(in)::b(:);real(dp),allocatable::x(:)
real(dp),allocatable::bm(:,:),sol(:,:),work(:)
if(size(b)/=f%n)error stop 'spam_solve: dimension mismatch'
allocate(bm(f%n,1),sol(f%n,1),work(f%n));bm(:,1)=b
call backsolves(f%n,f%nsuper,1,f%colindices,f%colpointers,f%entries,f%rowpointers,f%invpivot,f%pivot,f%supernodes,work,sol,bm)
allocate(x(f%n));x=sol(:,1)
end function spam_solve_vec

function spam_solve_mat(f,b) result(x)
type(spam_chol),intent(in)::f;real(dp),intent(in)::b(:,:);real(dp),allocatable::x(:,:)
integer::j
if(size(b,1)/=f%n)error stop 'spam_solve: dimension mismatch'
allocate(x(f%n,size(b,2)))
do j=1,size(b,2)
   x(:,j)=spam_solve_vec(f,b(:,j))
end do
end function spam_solve_mat

function spam_forwardsolve(f,b) result(x)
type(spam_chol),intent(in)::f;real(dp),intent(in)::b(:,:);real(dp),allocatable::x(:,:)
real(dp),allocatable::work(:),bm(:,:),sol(:,:);integer::j
if(size(b,1)/=f%n)error stop 'spam_forwardsolve: dimension mismatch'
allocate(x(f%n,size(b,2)),work(f%n),bm(f%n,1),sol(f%n,1))
do j=1,size(b,2)
   bm(:,1)=b(:,j)
   call pivotforwardsolve(f%n,f%nsuper,1,f%colindices,f%colpointers,f%entries,f%rowpointers, &
        f%invpivot,f%pivot,f%supernodes,work,sol,bm)
   x(:,j)=sol(:,1)
end do
end function spam_forwardsolve

function spam_backsolve(f,b) result(x)
type(spam_chol),intent(in)::f;real(dp),intent(in)::b(:,:);real(dp),allocatable::x(:,:)
real(dp),allocatable::work(:),bm(:,:),sol(:,:);integer::j
if(size(b,1)/=f%n)error stop 'spam_backsolve: dimension mismatch'
allocate(x(f%n,size(b,2)),work(f%n),bm(f%n,1),sol(f%n,1))
do j=1,size(b,2)
   bm(:,1)=b(:,j)
   call pivotbacksolve(f%n,f%nsuper,1,f%colindices,f%colpointers,f%entries,f%rowpointers, &
        f%invpivot,f%pivot,f%supernodes,work,sol,bm)
   x(:,j)=sol(:,1)
end do
end function spam_backsolve

real(dp) function spam_logdet(f) result(v)
type(spam_chol),intent(in)::f;integer::i
if(f%info/=0)then;v=-huge(1.0_dp);return;end if
v=0.0_dp
do i=1,f%n
 if(f%entries(f%rowpointers(i))<=0.0_dp)then;v=-huge(1.0_dp);return;end if
 v=v+2.0_dp*log(f%entries(f%rowpointers(i)))
end do
end function spam_logdet
real(dp) function spam_determinant(f) result(v)
type(spam_chol),intent(in)::f;v=exp(spam_logdet(f));end function spam_determinant
real(dp) function chol_factor_logdet(f) result(v)
type(spam_chol),intent(in)::f
v=0.5_dp*spam_logdet(f)
end function chol_factor_logdet

function spam_chol2inv(f,eps) result(a)
use spam_csr, only: csr_from_dense
type(spam_chol),intent(in)::f
real(dp),intent(in),optional::eps
type(csr_matrix)::a
real(dp),allocatable::eye(:,:),inv(:,:)
integer::i
allocate(eye(f%n,f%n));eye=0.0_dp
do i=1,f%n;eye(i,i)=1.0_dp;end do
inv=spam_solve_mat(f,eye)
a=csr_from_dense(inv,eps)
end function spam_chol2inv

function chol_to_dense(f) result(l)
type(spam_chol),intent(in)::f;real(dp),allocatable::l(:,:)
integer::s,j,base,k,lenrow,firstcol,offset,row
allocate(l(f%n,f%n));l=0.0_dp
! Ng-Peyton stores each supernode column densely in lnz; reconstruct the permuted lower factor.
do s=1,f%nsuper
   firstcol=f%supernodes(s); lenrow=f%colpointers(s+1)-f%colpointers(s)
   do j=firstcol,f%supernodes(s+1)-1
      base=f%rowpointers(j); offset=j-firstcol
      do k=offset,lenrow-1
         row=f%colindices(f%colpointers(s)+k)
         l(row,j)=f%entries(base+k-offset)
      end do
   end do
end do
end function chol_to_dense

pure function lower(s) result(t)
character(len=*),intent(in)::s;character(len=len(s))::t;integer::i,c
t=s;do i=1,len(s);c=iachar(t(i:i));if(c>=65.and.c<=90)t(i:i)=achar(c+32);end do
end function lower
end module spam_cholesky
