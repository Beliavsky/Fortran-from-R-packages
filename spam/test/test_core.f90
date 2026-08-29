program test_core
use spam
implicit none
integer,parameter::n=4
real(dp)::d(3,3),expect(3,3),expect23(2,3),v(5),h(3),th(3),loc(3,1),rw1(4,4),rw2(5,5)
integer::ir(5),jc(5)
type(csr_matrix)::a,c,q,dist
real(dp),allocatable::x(:,:),cv(:)
integer::i
ir=[1,1,1,2,8]
jc=[1,1,3,2,1]
v=[1d0,2d0,4d0,5d0,9d0]
a=csr_from_triplet(2,3,ir,jc,v)
call check(a%valid(),'triplet valid')
call check(a%nnz()==3,'duplicates/out of range')
x=csr_to_dense(a)
expect23=0
expect23(1,1)=3
expect23(1,3)=4
expect23(2,2)=5
call close2(x,expect23,1d-14,'triplet values')
c=csr_circulant([1d0,2d0,3d0])
x=csr_to_dense(c)
expect=reshape([1d0,3d0,2d0,2d0,1d0,3d0,3d0,2d0,1d0],[3,3])
call close2(x,expect,1d-14,'circulant')
q=precmat_rw1(4)
rw1=reshape([1d0,-1d0,0d0,0d0,-1d0,2d0,-1d0,0d0,0d0,-1d0,2d0,-1d0,0d0,0d0,-1d0,1d0],[4,4])
call close2(csr_to_dense(q),rw1,1d-14,'rw1')
q=precmat_rw2(5)
rw2=0
do i=1,3
 rw2(i:i+2,i:i+2)=rw2(i:i+2,i:i+2)+reshape([1d0,-2d0,1d0,-2d0,4d0,-2d0,1d0,-2d0,1d0],[3,3])
end do
call close2(csr_to_dense(q),rw2,1d-14,'rw2')
h=[0d0,0.5d0,1.5d0]
th=[1d0,2d0,0.25d0]
cv=cov_sph(h,th)
call close1(cv,[2.25d0,0.625d0,0d0],1d-14,'cov.sph')
loc(:,1)=[0d0,0.5d0,2d0]
dist=nearest_dist(loc,method='euclidean',delta=2d0,full=.true.)
call check(dist%nnz()==9,'nearest.dist keeps diagonal zeros')
call close2(csr_to_dense(cov_exp(dist,[1d0,2d0,0.5d0])), &
 reshape([2.5d0,2d0*exp(-0.5d0),2d0*exp(-2d0),2d0*exp(-0.5d0),2.5d0,2d0*exp(-1.5d0), &
          2d0*exp(-2d0),2d0*exp(-1.5d0),2.5d0],[3,3]),1d-13,'cov.exp')
print *,'test_core: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(*),intent(in)::msg
if(.not.ok)then
print *,msg
error stop
end if
end
subroutine close1(a,b,t,msg)
real(dp),intent(in)::a(:),b(:),t
character(*),intent(in)::msg
call check(maxval(abs(a-b))<=t,msg)
end
subroutine close2(a,b,t,msg)
real(dp),intent(in)::a(:,:),b(:,:),t
character(*),intent(in)::msg
call check(maxval(abs(a-b))<=t,msg)
end
end program
