module spam_block
use spam_kinds,only:dp
use spam_types,only:csr_matrix
implicit none
private
public::gmult
contains
function gmult(x,splits,fact) result(y)
type(csr_matrix),intent(in)::x
integer,intent(in)::splits(:)
real(dp),intent(in)::fact(:,:)
type(csr_matrix)::y
integer::i,k,rb,cb
if(size(fact,1)/=size(fact,2).or.size(fact,1)/=size(splits)-1)error stop 'gmult: factor/splits mismatch'
y=x
do i=1,x%nrow
 rb=block_index(i,splits)
 do k=x%rowpointers(i),x%rowpointers(i+1)-1
  cb=block_index(x%colindices(k),splits)
  y%entries(k)=x%entries(k)*fact(rb,cb)
 end do
end do
end function
integer function block_index(i,splits) result(b)
integer,intent(in)::i,splits(:);integer::k
b=size(splits)-1
do k=1,size(splits)-1
 if(i>=splits(k).and.i<splits(k+1))then;b=k;return;end if
end do
end function
end module spam_block
