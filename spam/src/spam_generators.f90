module spam_generators
use spam_kinds,only:dp
use spam_types,only:csr_matrix
use spam_csr,only:csr_from_triplet,csr_from_dense,csr_to_dense
use r_compat,only:runif_vec
implicit none
private
public::spam_random_matrix
abstract interface
 function random_values(n) result(x)
  import dp
  integer,intent(in)::n
  real(dp),allocatable::x(:)
 end function
end interface
contains
function spam_random_matrix(nrow,ncol,density,sym,spd,distribution) result(a)
integer,intent(in)::nrow
integer,intent(in),optional::ncol
real(dp),intent(in),optional::density
logical,intent(in),optional::sym,spd
procedure(random_values),optional::distribution
type(csr_matrix)::a
integer::nc,i,j,k,nz
real(dp)::rho
logical::sy,pd
real(dp),allocatable::u(:),v(:),d(:,:)
integer,allocatable::ir(:),jc(:)
nc=nrow
if(present(ncol))nc=ncol
rho=0.5_dp
if(present(density))rho=density
sy=.false.
if(present(sym))sy=sym
pd=.false.
if(present(spd))pd=spd
if(nrow<0.or.nc<0.or.rho<0.or.rho>1)error stop 'spam_random_matrix: invalid dimensions/density'
if((sy.or.pd).and.nrow/=nc)error stop 'spam_random_matrix: symmetric/SPD matrix must be square'
if(pd.and.rho==0)error stop 'spam_random_matrix: SPD requires positive density'
allocate(u(max(1,nrow*nc)))
u=runif_vec(size(u))
nz=count(u(:nrow*nc)<=rho)
if(nz==0)then
 allocate(d(nrow,nc))
 d=0.0_dp
 if(sy)then
 do i=1,min(nrow,nc)
 d(i,i)=1.0_dp
 end do
 end if
 a=csr_from_dense(d)
 return
end if
allocate(ir(nz),jc(nz),v(nz))
k=0
do i=1,nrow
do j=1,nc
 if(u((i-1)*nc+j)<=rho)then
 k=k+1
 ir(k)=i
 jc(k)=j
 v(k)=1.0_dp
 end if
end do
end do
if(present(distribution))v=distribution(nz)
a=csr_from_triplet(nrow,nc,ir,jc,v)
if(sy.or.pd)then
 d=csr_to_dense(a)
 do i=1,nrow
  do j=i+1,nrow
   if(d(i,j)/=0.0_dp)then
   d(j,i)=d(i,j)
   else
   d(i,j)=d(j,i)
   end if
  end do
 end do
 if(pd)then
  do i=1,nrow
   d(i,i)=d(i,i)+1.0_dp
  end do
  do i=1,nrow
   if(2.0_dp*abs(d(i,i))<=sum(abs(d(i,:))))then
    d(i,i)=sign(sum(abs(d(i,:)))+1.0_dp,merge(d(i,i),1.0_dp,d(i,i)/=0.0_dp))
   end if
  end do
 end if
 a=csr_from_dense(d)
end if
end function
end module spam_generators
