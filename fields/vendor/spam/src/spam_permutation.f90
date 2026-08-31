module spam_permutation
use spam_kinds,only:dp
use spam_types,only:csr_matrix
use spam_csr,only:csr_from_triplet
implicit none
private
public::permute_csr,inverse_permutation
contains
function inverse_permutation(p) result(inv)
integer,intent(in)::p(:);integer,allocatable::inv(:);integer::i
allocate(inv(size(p)));inv=0
do i=1,size(p)
 if(p(i)<1.or.p(i)>size(p).or.inv(p(i))/=0)error stop 'inverse_permutation: invalid permutation'
 inv(p(i))=i
end do
end function
function permute_csr(a,p,q,index_mode) result(b)
type(csr_matrix),intent(in)::a;integer,intent(in),optional::p(:),q(:);logical,intent(in),optional::index_mode
type(csr_matrix)::b
integer,allocatable::rp(:),cq(:),ir(:),jc(:);real(dp),allocatable::v(:);integer::i,k,nz,newr,newc
logical::ind
if(.not.present(p).and..not.present(q))error stop 'permute_csr: P or Q required'
ind=.false.;if(present(index_mode))ind=index_mode
if(present(p))then
 if(size(p)/=a%nrow)error stop 'permute_csr: P length mismatch'
 if(ind)then;rp=p;else;rp=inverse_permutation(p);end if
else
 allocate(rp(a%nrow));do i=1,a%nrow;rp(i)=i;end do
end if
if(present(q))then
 if(size(q)/=a%ncol)error stop 'permute_csr: Q length mismatch'
 if(ind)then;cq=q;else;cq=inverse_permutation(q);end if
else
 allocate(cq(a%ncol));do i=1,a%ncol;cq(i)=i;end do
end if
nz=a%nnz();allocate(ir(nz),jc(nz),v(nz));k=0
do i=1,a%nrow
 do newr=1,a%nrow
  if(rp(newr)==i)exit
 end do
 do nz=a%rowpointers(i),a%rowpointers(i+1)-1
  k=k+1;ir(k)=newr
  newc=0
  do newc=1,a%ncol;if(cq(newc)==a%colindices(nz))exit;end do
  jc(k)=newc;v(k)=a%entries(nz)
 end do
end do
b=csr_from_triplet(a%nrow,a%ncol,ir,jc,v,eps=-1.0_dp)
end function
end module spam_permutation
