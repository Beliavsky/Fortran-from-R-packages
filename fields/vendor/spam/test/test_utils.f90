program test_utils
use spam
implicit none
real(dp)::x(3,2),b(5,5),fact(2,2),dd
integer::p(3),splits(3)
type(csr_matrix)::a,c,q
real(dp),allocatable::z(:,:)
x=reshape([0d0,3d0,0d0,0d0,4d0,1d0],[3,2])
a=nearest_dist(x,method='euclidean',delta=10d0,full=.true.)
z=csr_to_dense(a);call check(abs(z(1,2)-5d0)<1d-14,'euclidean')
a=nearest_dist(x,method='maximum',delta=10d0,full=.true.)
z=csr_to_dense(a);call check(abs(z(1,2)-4d0)<1d-14,'maximum')
a=nearest_dist(x,method='minkowski',delta=10d0,p=1d0,full=.true.)
z=csr_to_dense(a);call check(abs(z(1,2)-7d0)<1d-14,'minkowski')
x=reshape([0d0,0d0,1d0,0d0,0d0,0d0],[3,2])
a=nearest_dist(x,method='greatcircle',delta=2d0,radius=6378.388d0,full=.true.)
z=csr_to_dense(a);dd=acos(-1d0)/180d0*6378.388d0
call check(abs(z(1,3)-dd)<1d-9,'greatcircle')
b=0; b(1,1)=1;b(2,2)=2;b(3,3)=3;b(4,4)=4;b(5,5)=5
a=csr_from_dense(b);p=[3,1,2]
! permutation test on a 3x3 principal block
a=csr_subset(a,[1,2,3],[1,2,3]);c=permute_csr(a,p=p,q=p,index_mode=.true.)
z=csr_to_dense(c);call check(all(abs(diagonal_dense(z)-[3d0,1d0,2d0])<1d-14),'permutation')
splits=[1,3,4];fact=reshape([2d0,4d0,3d0,5d0],[2,2]);c=gmult(a,splits,fact);z=csr_to_dense(c)
call check(abs(z(1,1)-2d0)<1d-14 .and. abs(z(3,3)-15d0)<1d-14,'gmult')
q=precmat_season(6,3);z=csr_to_dense(q);b=0;b(1:5,1:5)=0
call check(abs(z(1,1)-1d0)<1d-14 .and. abs(z(3,3)-3d0)<1d-14 .and. abs(z(1,3)-1d0)<1d-14,'season')
print *,'test_utils: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(*),intent(in)::msg
if(.not.ok)then;print *,msg;error stop;end if
end
function diagonal_dense(a) result(d)
real(dp),intent(in)::a(:,:);real(dp)::d(min(size(a,1),size(a,2)));integer::i
do i=1,size(d);d(i)=a(i,i);end do
end function
end program
