module spam_apply
use spam_kinds,only:dp
use spam_types,only:csr_matrix
use spam_csr,only:csr_from_triplet,csr_to_dense
implicit none
private
public::csr_map_entries,csr_apply_margin,csr_upper,csr_lower,csr_var
abstract interface
 function scalar_map(x) result(y)
  import dp
  real(dp),intent(in)::x
  real(dp)::y
 end function
 function vector_reduce(x) result(y)
  import dp
  real(dp),intent(in)::x(:)
  real(dp)::y
 end function
end interface
contains
function csr_map_entries(a,fun) result(b)
type(csr_matrix),intent(in)::a
procedure(scalar_map)::fun
type(csr_matrix)::b
integer::k
b=a
do k=1,b%nnz()
b%entries(k)=fun(b%entries(k))
if(.not.isfinite(b%entries(k)))b%entries(k)=0.0_dp
end do
end function
function csr_apply_margin(a,margin,fun) result(v)
type(csr_matrix),intent(in)::a
integer,intent(in)::margin
procedure(vector_reduce)::fun
real(dp),allocatable::v(:)
real(dp),allocatable::tmp(:)
integer::i,k,n
select case(margin)
case(1)
 allocate(v(a%nrow))
 do i=1,a%nrow
  n=a%rowpointers(i+1)-a%rowpointers(i)
  allocate(tmp(n))
  if(n>0)tmp=a%entries(a%rowpointers(i):a%rowpointers(i+1)-1)
  v(i)=fun(tmp)
  deallocate(tmp)
 end do
case(2)
 allocate(v(a%ncol))
 do i=1,a%ncol
 n=count(a%colindices==i)
 allocate(tmp(n))
 k=0
  if(n>0)then
  do k=1,a%nnz()
  if(a%colindices(k)==i)tmp(count(a%colindices(:k)==i))=a%entries(k)
  end do
  end if
  v(i)=fun(tmp)
  deallocate(tmp)
 end do
case default;error stop 'csr_apply_margin: margin must be 1 or 2'
end select
end function
function csr_upper(a,diag) result(b)
type(csr_matrix),intent(in)::a
logical,intent(in),optional::diag
type(csr_matrix)::b
b=tri_part(a,.true.,diag)
end function
function csr_lower(a,diag) result(b)
type(csr_matrix),intent(in)::a
logical,intent(in),optional::diag
type(csr_matrix)::b
b=tri_part(a,.false.,diag)
end function
function tri_part(a,upper,diag) result(b)
type(csr_matrix),intent(in)::a
logical,intent(in)::upper
logical,intent(in),optional::diag
type(csr_matrix)::b
integer,allocatable::ir(:),jc(:)
real(dp),allocatable::v(:)
integer::i,k,nz
logical::dg,keep
dg=.false.
if(present(diag))dg=diag
allocate(ir(a%nnz()),jc(a%nnz()),v(a%nnz()))
nz=0
do i=1,a%nrow
do k=a%rowpointers(i),a%rowpointers(i+1)-1
 if(upper)then
  keep=a%colindices(k)>i.or.(dg.and.a%colindices(k)==i)
 else
  keep=a%colindices(k)<i.or.(dg.and.a%colindices(k)==i)
 end if
 if(keep)then
 nz=nz+1
 ir(nz)=i
 jc(nz)=a%colindices(k)
 v(nz)=a%entries(k)
 end if
end do
end do
b=csr_from_triplet(a%nrow,a%ncol,ir(:nz),jc(:nz),v(:nz),eps=-1.0_dp)
end function
function csr_var(a) result(v)
type(csr_matrix),intent(in)::a
real(dp),allocatable::v(:,:)
real(dp),allocatable::x(:,:),mu(:)
integer::n,i
x=csr_to_dense(a)
n=a%nrow
if(n<2)error stop 'csr_var: at least two rows required'
allocate(mu(a%ncol))
mu=sum(x,dim=1)/real(n,dp)
allocate(v(a%ncol,a%ncol))
v=0.0_dp
do i=1,n
v=v+outer(x(i,:)-mu,x(i,:)-mu)
end do
v=v/real(n-1,dp)
end function
pure function outer(a,b) result(c)
real(dp),intent(in)::a(:),b(:)
real(dp)::c(size(a),size(b))
integer::i
do i=1,size(a)
c(i,:)=a(i)*b
end do
end function
pure logical function isfinite(x) result(ok)
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
real(dp),intent(in)::x
ok=ieee_is_finite(x)
end function
end module spam_apply
