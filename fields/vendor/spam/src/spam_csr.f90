module spam_csr
use spam_kinds, only: dp
use spam_types, only: csr_matrix
implicit none
private
public :: csr_from_dense, csr_from_triplet, csr_to_dense, csr_zeros, csr_diag, csr_identity
public :: csr_transpose, csr_add, csr_subtract, csr_scale, csr_hadamard, csr_matmul
public :: csr_matvec, csr_matmat, csr_kronecker, csr_rbind, csr_cbind, csr_subset
public :: csr_get, csr_diagonal, csr_set_diagonal, csr_row_sums, csr_col_sums
public :: csr_row_means, csr_col_means, csr_norm, csr_clean, csr_is_symmetric
public :: csr_diff, csr_toeplitz, csr_circulant, csr_crossprod, csr_tcrossprod
public :: csr_bdiag, csr_bandwidth

contains

function csr_zeros(nrow,ncol) result(a)
integer, intent(in) :: nrow,ncol
type(csr_matrix) :: a
a%nrow=nrow; a%ncol=ncol
allocate(a%entries(0),a%colindices(0),a%rowpointers(nrow+1))
a%rowpointers=1
end function csr_zeros

function csr_from_dense(x,eps) result(a)
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: eps
type(csr_matrix) :: a
real(dp) :: tol
integer :: i,j,k,nz
tol=0.0_dp; if(present(eps)) tol=eps
nz=count(abs(x)>tol)
a%nrow=size(x,1); a%ncol=size(x,2)
allocate(a%entries(nz),a%colindices(nz),a%rowpointers(a%nrow+1))
k=0; a%rowpointers(1)=1
do i=1,a%nrow
   do j=1,a%ncol
      if(abs(x(i,j))>tol) then
         k=k+1; a%entries(k)=x(i,j); a%colindices(k)=j
      end if
   end do
   a%rowpointers(i+1)=k+1
end do
end function csr_from_dense

function csr_from_triplet(nrow,ncol,ir,jc,val,eps) result(a)
integer, intent(in) :: nrow,ncol,ir(:),jc(:)
real(dp), intent(in) :: val(:)
real(dp), intent(in), optional :: eps
type(csr_matrix) :: a
integer :: i,k,p,nz,cap,c,nr,nvalid
integer, allocatable :: countrow(:),nextp(:),ord(:),cols(:),mark(:),rowcols(:)
real(dp), allocatable :: vals(:),work(:)
real(dp) :: tol
if(size(ir)/=size(jc) .or. size(ir)/=size(val)) error stop 'csr_from_triplet: size mismatch'
tol=0.0_dp; if(present(eps)) tol=eps
if(nrow<0 .or. ncol<0) error stop 'csr_from_triplet: negative dimension'
if(nrow==0 .or. ncol==0) then
   a=csr_zeros(nrow,ncol); return
end if
allocate(countrow(nrow),nextp(nrow))
countrow=0
! Ignore out-of-range triplets, matching upstream triplet3csr behavior.
do k=1,size(ir)
   if(ir(k)>=1 .and. ir(k)<=nrow .and. jc(k)>=1 .and. jc(k)<=ncol) countrow(ir(k))=countrow(ir(k))+1
end do
nvalid=sum(countrow)
allocate(ord(nvalid))
nextp(1)=1
do i=2,nrow
   nextp(i)=nextp(i-1)+countrow(i-1)
end do
countrow=nextp
do k=1,size(ir)
   if(ir(k)>=1 .and. ir(k)<=nrow .and. jc(k)>=1 .and. jc(k)<=ncol) then
      ord(countrow(ir(k)))=k; countrow(ir(k))=countrow(ir(k))+1
   end if
end do
cap=max(16,min(max(1,size(val)),1024))
allocate(vals(cap),cols(cap),mark(ncol),work(ncol),rowcols(max(1,ncol)))
mark=0; work=0.0_dp; nz=0
a%nrow=nrow; a%ncol=ncol
allocate(a%rowpointers(nrow+1)); a%rowpointers(1)=1
p=1
do i=1,nrow
   nr=0
   do while(p<=nvalid)
      k=ord(p)
      if(k<=0) then
         p=p+1; cycle
      end if
      if(ir(k)/=i) exit
      c=jc(k)
      if(mark(c)/=i) then
         nr=nr+1; rowcols(nr)=c; mark(c)=i; work(c)=val(k)
      else
         work(c)=work(c)+val(k)
      end if
      p=p+1
   end do
   call sort_int(rowcols,nr)
   do k=1,nr
      c=rowcols(k)
      if(abs(work(c))>tol) then
         nz=nz+1
         call ensure_capacity(vals,cols,nz)
         vals(nz)=work(c); cols(nz)=c
      end if
      work(c)=0.0_dp
   end do
   a%rowpointers(i+1)=nz+1
end do
allocate(a%entries(nz),a%colindices(nz))
if(nz>0) then
   a%entries=vals(:nz); a%colindices=cols(:nz)
end if
end function csr_from_triplet

function csr_to_dense(a) result(x)
type(csr_matrix), intent(in) :: a
real(dp), allocatable :: x(:,:)
integer :: i,k
allocate(x(a%nrow,a%ncol)); x=0.0_dp
do i=1,a%nrow
   do k=a%rowpointers(i),a%rowpointers(i+1)-1
      x(i,a%colindices(k))=a%entries(k)
   end do
end do
end function csr_to_dense

function csr_diag(d,nrow,ncol,eps) result(a)
real(dp), intent(in) :: d(:)
integer, intent(in), optional :: nrow,ncol
real(dp), intent(in), optional :: eps
type(csr_matrix) :: a
integer :: nr,nc,n,i,k
real(dp)::tol
nr=size(d); if(present(nrow)) nr=nrow
nc=nr; if(present(ncol)) nc=ncol
n=min(nr,nc); tol=0.0_dp; if(present(eps)) tol=eps
a%nrow=nr; a%ncol=nc
allocate(a%rowpointers(nr+1)); a%rowpointers(1)=1
k=0
do i=1,nr
   if(i<=n .and. i<=size(d)) then
      if(abs(d(i))>tol) k=k+1
   end if
   a%rowpointers(i+1)=k+1
end do
allocate(a%entries(k),a%colindices(k)); k=0
do i=1,min(n,size(d))
   if(abs(d(i))>tol) then
      k=k+1; a%entries(k)=d(i); a%colindices(k)=i
   end if
end do
end function csr_diag

function csr_identity(n,scale) result(a)
integer,intent(in)::n
real(dp),intent(in),optional::scale
type(csr_matrix)::a
real(dp),allocatable::d(:)
real(dp)::s
s=1.0_dp; if(present(scale))s=scale
allocate(d(n));d=s
a=csr_diag(d)
end function csr_identity

function csr_transpose(a) result(at)
type(csr_matrix),intent(in)::a
type(csr_matrix)::at
integer::i,k,j,nz
integer,allocatable::cnt(:),next(:)
nz=a%nnz(); at%nrow=a%ncol; at%ncol=a%nrow
allocate(at%entries(nz),at%colindices(nz),at%rowpointers(at%nrow+1),cnt(at%nrow),next(at%nrow))
cnt=0
do k=1,nz; cnt(a%colindices(k))=cnt(a%colindices(k))+1; end do
at%rowpointers(1)=1
do i=1,at%nrow; at%rowpointers(i+1)=at%rowpointers(i)+cnt(i); end do
next=at%rowpointers(:at%nrow)
do i=1,a%nrow
   do k=a%rowpointers(i),a%rowpointers(i+1)-1
      j=a%colindices(k); at%entries(next(j))=a%entries(k); at%colindices(next(j))=i; next(j)=next(j)+1
   end do
end do
end function csr_transpose

function csr_add(a,b,eps) result(c)
type(csr_matrix),intent(in)::a,b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
real(dp)::tol,v
integer::i,ka,kb,eaa,ebb,nz,cap,col
real(dp),allocatable::vals(:)
integer,allocatable::cols(:)
if(a%nrow/=b%nrow .or. a%ncol/=b%ncol) error stop 'csr_add: nonconformable'
tol=0.0_dp;if(present(eps))tol=eps
cap=max(16,a%nnz()+b%nnz());allocate(vals(cap),cols(cap),c%rowpointers(a%nrow+1))
c%nrow=a%nrow;c%ncol=a%ncol;nz=0;c%rowpointers(1)=1
do i=1,a%nrow
   ka=a%rowpointers(i); eaa=a%rowpointers(i+1)-1
   kb=b%rowpointers(i); ebb=b%rowpointers(i+1)-1
   do while(ka<=eaa .or. kb<=ebb)
      if(kb>ebb .or. (ka<=eaa .and. a%colindices(ka)<b%colindices(kb))) then
         col=a%colindices(ka);v=a%entries(ka);ka=ka+1
      else if(ka>eaa .or. b%colindices(kb)<a%colindices(ka)) then
         col=b%colindices(kb);v=b%entries(kb);kb=kb+1
      else
         col=a%colindices(ka);v=a%entries(ka)+b%entries(kb);ka=ka+1;kb=kb+1
      end if
      if(abs(v)>tol) then
         nz=nz+1;call ensure_capacity(vals,cols,nz);vals(nz)=v;cols(nz)=col
      end if
   end do
   c%rowpointers(i+1)=nz+1
end do
allocate(c%entries(nz),c%colindices(nz));if(nz>0)then;c%entries=vals(:nz);c%colindices=cols(:nz);end if
end function csr_add

function csr_subtract(a,b,eps) result(c)
type(csr_matrix),intent(in)::a,b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
type(csr_matrix)::mb
mb=csr_scale(b,-1.0_dp);c=csr_add(a,mb,eps)
end function csr_subtract

function csr_scale(a,s,eps) result(c)
type(csr_matrix),intent(in)::a
real(dp),intent(in)::s
real(dp),intent(in),optional::eps
type(csr_matrix)::c
real(dp)::tol
integer::i,k,nz
tol=0.0_dp;if(present(eps))tol=eps
c%nrow=a%nrow;c%ncol=a%ncol
nz=count(abs(s*a%entries)>tol)
allocate(c%entries(nz),c%colindices(nz),c%rowpointers(c%nrow+1));k=0;c%rowpointers(1)=1
do i=1,a%nrow
   do nz=a%rowpointers(i),a%rowpointers(i+1)-1
      if(abs(s*a%entries(nz))>tol)then;k=k+1;c%entries(k)=s*a%entries(nz);c%colindices(k)=a%colindices(nz);end if
   end do
   c%rowpointers(i+1)=k+1
end do
end function csr_scale

function csr_hadamard(a,b,eps) result(c)
type(csr_matrix),intent(in)::a,b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
real(dp)::tol,v
integer::i,ka,kb,nz,cap
real(dp),allocatable::vals(:);integer,allocatable::cols(:)
if(a%nrow/=b%nrow .or. a%ncol/=b%ncol)error stop 'csr_hadamard: nonconformable'
tol=0.0_dp;if(present(eps))tol=eps
cap=max(16,min(a%nnz(),b%nnz()));allocate(vals(cap),cols(cap),c%rowpointers(a%nrow+1));nz=0
c%nrow=a%nrow;c%ncol=a%ncol;c%rowpointers(1)=1
do i=1,a%nrow
   ka=a%rowpointers(i);kb=b%rowpointers(i)
   do while(ka<a%rowpointers(i+1) .and. kb<b%rowpointers(i+1))
      if(a%colindices(ka)<b%colindices(kb))then;ka=ka+1
      else if(b%colindices(kb)<a%colindices(ka))then;kb=kb+1
      else
         v=a%entries(ka)*b%entries(kb)
         if(abs(v)>tol)then;nz=nz+1;call ensure_capacity(vals,cols,nz);vals(nz)=v;cols(nz)=a%colindices(ka);end if
         ka=ka+1;kb=kb+1
      end if
   end do
   c%rowpointers(i+1)=nz+1
end do
allocate(c%entries(nz),c%colindices(nz));if(nz>0)then;c%entries=vals(:nz);c%colindices=cols(:nz);end if
end function csr_hadamard

function csr_matmul(a,b,eps) result(c)
type(csr_matrix),intent(in)::a,b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
real(dp)::tol,v
integer::i,ka,kb,j,col,nr,nz,cap
integer,allocatable::mark(:),rowcols(:),cols(:)
real(dp),allocatable::work(:),vals(:)
if(a%ncol/=b%nrow)error stop 'csr_matmul: nonconformable'
tol=0.0_dp;if(present(eps))tol=eps
allocate(mark(b%ncol),work(b%ncol),rowcols(max(1,b%ncol)));mark=0;work=0.0_dp
cap=max(16,min(max(1,a%nnz()+b%nnz()),max(16,a%nrow*4)));allocate(vals(cap),cols(cap),c%rowpointers(a%nrow+1))
c%nrow=a%nrow;c%ncol=b%ncol;nz=0;c%rowpointers(1)=1
do i=1,a%nrow
   nr=0
   do ka=a%rowpointers(i),a%rowpointers(i+1)-1
      j=a%colindices(ka)
      do kb=b%rowpointers(j),b%rowpointers(j+1)-1
         col=b%colindices(kb)
         if(mark(col)/=i)then;nr=nr+1;rowcols(nr)=col;mark(col)=i;work(col)=a%entries(ka)*b%entries(kb)
         else;work(col)=work(col)+a%entries(ka)*b%entries(kb);end if
      end do
   end do
   call sort_int(rowcols,nr)
   do j=1,nr
      col=rowcols(j);v=work(col)
      if(abs(v)>tol)then;nz=nz+1;call ensure_capacity(vals,cols,nz);vals(nz)=v;cols(nz)=col;end if
      work(col)=0.0_dp
   end do
   c%rowpointers(i+1)=nz+1
end do
allocate(c%entries(nz),c%colindices(nz));if(nz>0)then;c%entries=vals(:nz);c%colindices=cols(:nz);end if
end function csr_matmul

function csr_matvec(a,x) result(y)
type(csr_matrix),intent(in)::a
real(dp),intent(in)::x(:)
real(dp),allocatable::y(:)
integer::i,k
if(size(x)/=a%ncol)error stop 'csr_matvec: nonconformable'
allocate(y(a%nrow));y=0.0_dp
do i=1,a%nrow
   do k=a%rowpointers(i),a%rowpointers(i+1)-1;y(i)=y(i)+a%entries(k)*x(a%colindices(k));end do
end do
end function csr_matvec

function csr_matmat(a,x) result(y)
type(csr_matrix),intent(in)::a
real(dp),intent(in)::x(:,:)
real(dp),allocatable::y(:,:)
integer::i,k,j
if(size(x,1)/=a%ncol)error stop 'csr_matmat: nonconformable'
allocate(y(a%nrow,size(x,2)));y=0.0_dp
do i=1,a%nrow
 do k=a%rowpointers(i),a%rowpointers(i+1)-1
   do j=1,size(x,2);y(i,j)=y(i,j)+a%entries(k)*x(a%colindices(k),j);end do
 end do
end do
end function csr_matmat

function csr_kronecker(a,b,eps) result(c)
type(csr_matrix),intent(in)::a,b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
real(dp)::tol,v
integer::ia,ka,ib,kb,row,nz,cap,col
real(dp),allocatable::vals(:);integer,allocatable::cols(:)
tol=0.0_dp;if(present(eps))tol=eps
cap=max(16,a%nnz()*max(1,b%nnz()));allocate(vals(cap),cols(cap),c%rowpointers(a%nrow*b%nrow+1))
c%nrow=a%nrow*b%nrow;c%ncol=a%ncol*b%ncol;nz=0;c%rowpointers(1)=1;row=0
do ia=1,a%nrow
 do ib=1,b%nrow
   row=row+1
   do ka=a%rowpointers(ia),a%rowpointers(ia+1)-1
    do kb=b%rowpointers(ib),b%rowpointers(ib+1)-1
      v=a%entries(ka)*b%entries(kb)
      if(abs(v)>tol)then
        nz=nz+1;call ensure_capacity(vals,cols,nz);vals(nz)=v
        col=(a%colindices(ka)-1)*b%ncol+b%colindices(kb);cols(nz)=col
      end if
    end do
   end do
   c%rowpointers(row+1)=nz+1
 end do
end do
allocate(c%entries(nz),c%colindices(nz));if(nz>0)then;c%entries=vals(:nz);c%colindices=cols(:nz);end if
end function csr_kronecker

function csr_rbind(a,b) result(c)
type(csr_matrix),intent(in)::a,b
type(csr_matrix)::c
integer::nza,nzb
if(a%ncol/=b%ncol)error stop 'csr_rbind: ncol mismatch'
nza=a%nnz();nzb=b%nnz();c%nrow=a%nrow+b%nrow;c%ncol=a%ncol
allocate(c%entries(nza+nzb),c%colindices(nza+nzb),c%rowpointers(c%nrow+1))
if(nza>0)then;c%entries(:nza)=a%entries;c%colindices(:nza)=a%colindices;end if
if(nzb>0)then;c%entries(nza+1:)=b%entries;c%colindices(nza+1:)=b%colindices;end if
c%rowpointers(:a%nrow+1)=a%rowpointers
c%rowpointers(a%nrow+2:)=b%rowpointers(2:)+nza
end function csr_rbind

function csr_cbind(a,b) result(c)
type(csr_matrix),intent(in)::a,b
type(csr_matrix)::c
integer::i,ka,kb,nz
if(a%nrow/=b%nrow)error stop 'csr_cbind: nrow mismatch'
c%nrow=a%nrow;c%ncol=a%ncol+b%ncol;nz=a%nnz()+b%nnz()
allocate(c%entries(nz),c%colindices(nz),c%rowpointers(c%nrow+1));nz=0;c%rowpointers(1)=1
do i=1,c%nrow
 do ka=a%rowpointers(i),a%rowpointers(i+1)-1;nz=nz+1;c%entries(nz)=a%entries(ka);c%colindices(nz)=a%colindices(ka);end do
 do kb=b%rowpointers(i),b%rowpointers(i+1)-1;nz=nz+1;c%entries(nz)=b%entries(kb);c%colindices(nz)=b%colindices(kb)+a%ncol;end do
 c%rowpointers(i+1)=nz+1
end do
end function csr_cbind

function csr_subset(a,rows,cols,eps) result(c)
type(csr_matrix),intent(in)::a
integer,intent(in)::rows(:),cols(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::c
integer::i,j,k,nz,cap,cc
integer,allocatable::map(:),outc(:)
real(dp),allocatable::outv(:)
real(dp)::tol
tol=0.0_dp;if(present(eps))tol=eps
allocate(map(a%ncol));map=0
do j=1,size(cols);if(cols(j)>=1.and.cols(j)<=a%ncol)map(cols(j))=j;end do
cap=max(16,min(a%nnz(),max(1,size(rows)*4)));allocate(outv(cap),outc(cap),c%rowpointers(size(rows)+1));nz=0
c%nrow=size(rows);c%ncol=size(cols);c%rowpointers(1)=1
do i=1,size(rows)
 if(rows(i)>=1.and.rows(i)<=a%nrow)then
  do k=a%rowpointers(rows(i)),a%rowpointers(rows(i)+1)-1
   cc = map(a%colindices(k))
   if (cc > 0 .and. abs(a%entries(k)) > tol) then
      nz = nz + 1
      call ensure_capacity(outv,outc,nz)
      outv(nz) = a%entries(k)
      outc(nz) = cc
   end if
  end do
 end if
 c%rowpointers(i+1)=nz+1
end do
allocate(c%entries(nz),c%colindices(nz));if(nz>0)then;c%entries=outv(:nz);c%colindices=outc(:nz);end if
end function csr_subset

real(dp) function csr_get(a,i,j) result(v)
type(csr_matrix),intent(in)::a
integer,intent(in)::i,j
integer::k
v=0.0_dp;if(i<1.or.i>a%nrow.or.j<1.or.j>a%ncol)return
do k=a%rowpointers(i),a%rowpointers(i+1)-1
 if(a%colindices(k)==j)then;v=a%entries(k);return;else if(a%colindices(k)>j)then;return;end if
end do
end function csr_get

function csr_diagonal(a) result(d)
type(csr_matrix),intent(in)::a
real(dp),allocatable::d(:)
integer::i
allocate(d(min(a%nrow,a%ncol)))
do i=1,size(d);d(i)=csr_get(a,i,i);end do
end function csr_diagonal

function csr_set_diagonal(a,d,add,eps) result(c)
type(csr_matrix),intent(in)::a
real(dp),intent(in)::d(:)
logical,intent(in),optional::add
real(dp),intent(in),optional::eps
type(csr_matrix)::c
type(csr_matrix)::dm
real(dp),allocatable::dd(:)
logical::doadd
integer::n,i
n=min(a%nrow,a%ncol);allocate(dd(n));dd=0.0_dp;dd(:min(n,size(d)))=d(:min(n,size(d)))
doadd=.false.;if(present(add))doadd=add
if(doadd)then
 dm=csr_diag(dd,a%nrow,a%ncol,eps);c=csr_add(a,dm,eps)
else
 c=a
 ! easiest robust replacement: add (new-old) diagonal
 do i=1,n;dd(i)=dd(i)-csr_get(a,i,i);end do
 dm=csr_diag(dd,a%nrow,a%ncol,eps);c=csr_add(a,dm,eps)
end if
end function csr_set_diagonal

function csr_row_sums(a) result(s)
type(csr_matrix),intent(in)::a
real(dp),allocatable::s(:)
integer::i
allocate(s(a%nrow));do i=1,a%nrow;s(i)=sum(a%entries(a%rowpointers(i):a%rowpointers(i+1)-1));end do
end function csr_row_sums

function csr_col_sums(a) result(s)
type(csr_matrix),intent(in)::a
real(dp),allocatable::s(:)
integer::k
allocate(s(a%ncol));s=0.0_dp;do k=1,a%nnz();s(a%colindices(k))=s(a%colindices(k))+a%entries(k);end do
end function csr_col_sums

function csr_row_means(a,structure_based) result(s)
type(csr_matrix),intent(in)::a
logical,intent(in),optional::structure_based
real(dp),allocatable::s(:)
logical::sb;integer::i,n
sb=.false.;if(present(structure_based))sb=structure_based
s=csr_row_sums(a)
do i=1,a%nrow
 n=merge(a%rowpointers(i+1)-a%rowpointers(i),a%ncol,sb)
 if(n>0)s(i)=s(i)/real(n,dp)
end do
end function csr_row_means

function csr_col_means(a,structure_based) result(s)
type(csr_matrix),intent(in)::a
logical,intent(in),optional::structure_based
real(dp),allocatable::s(:)
integer,allocatable::cnt(:);logical::sb;integer::k,j
sb=.false.;if(present(structure_based))sb=structure_based
s=csr_col_sums(a);allocate(cnt(a%ncol));cnt=0
if(sb)then
 do k=1,a%nnz();cnt(a%colindices(k))=cnt(a%colindices(k))+1;end do
 do j=1,a%ncol;if(cnt(j)>0)s(j)=s(j)/real(cnt(j),dp);end do
else
 if(a%nrow>0)s=s/real(a%nrow,dp)
end if
end function csr_col_means

real(dp) function csr_norm(a,type) result(v)
type(csr_matrix),intent(in)::a
character(len=*),intent(in),optional::type
character(len=1)::t
real(dp),allocatable::s(:)
t='m';if(present(type))t=type(1:1)
select case(t)
case('o','O','1');s=csr_col_sums_abs(a);v=maxval(s)
case('i','I');s=csr_row_sums_abs(a);v=maxval(s)
case('f','F','h','H');v=sqrt(sum(a%entries*a%entries))
case default;if(a%nnz()>0)then;v=maxval(abs(a%entries));else;v=0.0_dp;end if
end select
end function csr_norm

function csr_clean(a,eps) result(c)
type(csr_matrix),intent(in)::a
real(dp),intent(in)::eps
type(csr_matrix)::c
integer::i,k,nz
c%nrow=a%nrow;c%ncol=a%ncol;nz=count(abs(a%entries)>eps)
allocate(c%entries(nz),c%colindices(nz),c%rowpointers(c%nrow+1));nz=0;c%rowpointers(1)=1
do i=1,a%nrow
 do k=a%rowpointers(i),a%rowpointers(i+1)-1
  if(abs(a%entries(k))>eps)then;nz=nz+1;c%entries(nz)=a%entries(k);c%colindices(nz)=a%colindices(k);end if
 end do
 c%rowpointers(i+1)=nz+1
end do
end function csr_clean

logical function csr_is_symmetric(a,tol) result(ok)
type(csr_matrix),intent(in)::a
real(dp),intent(in),optional::tol
real(dp)::t
integer::i,k
t=100.0_dp*epsilon(1.0_dp);if(present(tol))t=tol
if(a%nrow/=a%ncol)then;ok=.false.;return;end if
ok=.true.
do i=1,a%nrow
 do k=a%rowpointers(i),a%rowpointers(i+1)-1
  if(abs(a%entries(k)-csr_get(a,a%colindices(k),i))>t)then;ok=.false.;return;end if
 end do
end do
end function csr_is_symmetric

function csr_diff(a,lag,differences,eps) result(c)
type(csr_matrix),intent(in)::a
integer,intent(in),optional::lag,differences
real(dp),intent(in),optional::eps
type(csr_matrix)::c,d
integer::lg,nd,i
integer,allocatable::rows1(:),rows2(:),cols(:)
lg=1;if(present(lag))lg=lag;nd=1;if(present(differences))nd=differences
c=a
allocate(cols(a%ncol));cols=[(i,i=1,a%ncol)]
do i=1,nd
 if(c%nrow<=lg)then;c=csr_zeros(0,c%ncol);exit;end if
 allocate(rows1(c%nrow-lg),rows2(c%nrow-lg));rows1=[(i,i=1,c%nrow-lg)];rows2=rows1+lg
 d=csr_subtract(csr_subset(c,rows2,cols),csr_subset(c,rows1,cols),eps)
 call move_csr(d,c);deallocate(rows1,rows2)
end do
end function csr_diff

function csr_toeplitz(x,y,eps) result(a)
real(dp),intent(in)::x(:)
real(dp),intent(in),optional::y(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::a
real(dp),allocatable::d(:,:)
integer::n,i,j
n=size(x);if(present(y))then;if(size(y)/=n)error stop 'csr_toeplitz: y size';end if
allocate(d(n,n))
do i=1,n
 do j=1,n
  if(j>=i)then;d(i,j)=x(j-i+1)
  else;if(present(y))then;d(i,j)=y(i-j+1);else;d(i,j)=x(i-j+1);end if;end if
 end do
end do
a=csr_from_dense(d,eps)
end function csr_toeplitz

function csr_circulant(x,eps) result(a)
real(dp),intent(in)::x(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::a
real(dp),allocatable::d(:,:)
integer::n,i,j,idx
n=size(x);allocate(d(n,n))
do i=1,n;do j=1,n;idx=modulo(j-i,n)+1;d(i,j)=x(idx);end do;end do
a=csr_from_dense(d,eps)
end function csr_circulant

function csr_crossprod(a,b,eps) result(c)
type(csr_matrix),intent(in)::a
type(csr_matrix),intent(in),optional::b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
if(present(b))then;c=csr_matmul(csr_transpose(a),b,eps);else;c=csr_matmul(csr_transpose(a),a,eps);end if
end function csr_crossprod

function csr_tcrossprod(a,b,eps) result(c)
type(csr_matrix),intent(in)::a
type(csr_matrix),intent(in),optional::b
real(dp),intent(in),optional::eps
type(csr_matrix)::c
if(present(b))then;c=csr_matmul(a,csr_transpose(b),eps);else;c=csr_matmul(a,csr_transpose(a),eps);end if
end function csr_tcrossprod

function csr_bdiag(blocks) result(c)
type(csr_matrix),intent(in)::blocks(:)
type(csr_matrix)::c
integer::br,bc,b,k,i,nz,total
br=0;bc=0;total=0
do b=1,size(blocks);br=br+blocks(b)%nrow;bc=bc+blocks(b)%ncol;total=total+blocks(b)%nnz();end do
c%nrow=br;c%ncol=bc;allocate(c%entries(total),c%colindices(total),c%rowpointers(br+1));nz=0;br=0;bc=0;c%rowpointers(1)=1
do b=1,size(blocks)
 do i=1,blocks(b)%nrow
  do k=blocks(b)%rowpointers(i),blocks(b)%rowpointers(i+1)-1
   nz=nz+1;c%entries(nz)=blocks(b)%entries(k);c%colindices(nz)=blocks(b)%colindices(k)+bc
  end do
  br=br+1;c%rowpointers(br+1)=nz+1
 end do
 bc=bc+blocks(b)%ncol
end do
end function csr_bdiag

function csr_row_sums_abs(a) result(s)
type(csr_matrix),intent(in)::a
real(dp),allocatable::s(:);integer::i
allocate(s(a%nrow));do i=1,a%nrow;s(i)=sum(abs(a%entries(a%rowpointers(i):a%rowpointers(i+1)-1)));end do
end function csr_row_sums_abs
function csr_col_sums_abs(a) result(s)
type(csr_matrix),intent(in)::a
real(dp),allocatable::s(:);integer::k
allocate(s(a%ncol));s=0.0_dp;do k=1,a%nnz();s(a%colindices(k))=s(a%colindices(k))+abs(a%entries(k));end do
end function csr_col_sums_abs

subroutine ensure_capacity(vals,cols,need)
real(dp),allocatable,intent(inout)::vals(:)
integer,allocatable,intent(inout)::cols(:)
integer,intent(in)::need
real(dp),allocatable::v2(:);integer,allocatable::c2(:);integer::nc
if(need<=size(vals))return
nc=max(need,max(16,2*size(vals)));allocate(v2(nc),c2(nc));if(size(vals)>0)then;v2(:size(vals))=vals;c2(:size(cols))=cols;end if
call move_alloc(v2,vals);call move_alloc(c2,cols)
end subroutine ensure_capacity

subroutine sort_int(x,n)
integer,intent(inout)::x(:);integer,intent(in)::n
integer::i,j,t
! insertion sort is efficient for typical sparse rows and preserves determinism.
do i=2,n;t=x(i);j=i-1;do while(j>=1);if(x(j)<=t)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=t;end do
end subroutine sort_int

subroutine move_csr(src,dst)
type(csr_matrix),intent(inout)::src,dst
dst%nrow=src%nrow;dst%ncol=src%ncol
call move_alloc(src%entries,dst%entries)
call move_alloc(src%colindices,dst%colindices)
call move_alloc(src%rowpointers,dst%rowpointers)
end subroutine move_csr

function csr_bandwidth(a) result(bw)
type(csr_matrix),intent(in)::a
integer::bw(2),i,k,d
bw=0
do i=1,a%nrow
   do k=a%rowpointers(i),a%rowpointers(i+1)-1
      d=a%colindices(k)-i
      if(d<0)bw(1)=max(bw(1),-d)
      if(d>0)bw(2)=max(bw(2), d)
   end do
end do
end function csr_bandwidth

end module spam_csr
