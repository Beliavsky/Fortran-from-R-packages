program test_cholesky_eigen
use spam
implicit none
real(dp)::d(4,4),b(4)
real(dp),allocatable::x(:)
type(csr_matrix)::a
type(spam_chol)::f
type(eigen_result)::e
d=0
d(1,1)=4
d(2,2)=3
d(3,3)=2
d(4,4)=5
d(1,2)=1
d(2,1)=1
d(2,3)=0.5d0
d(3,2)=0.5d0
a=csr_from_dense(d)
f=spam_chol_factor(a,'mmd')
call check(f%info==0,'chol info')
b=[1d0,2d0,3d0,4d0]
x=spam_solve(f,b)
call check(maxval(abs(matmul(d,x)-b))<1d-11,'solve')
call check(abs(spam_logdet(f)-log(determinant4(d)))<1d-11,'logdet')
e=spam_eigen_symmetric(a,2,'LA',ncv=4)
call check(e%nconv==2,'eigen convergence')
call check(maxval(abs(matmul(d,e%vectors)-e%vectors*spread(e%values,1,4)))<1d-9,'eigen residual')
print *,'test_cholesky_eigen: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(*),intent(in)::msg
if(.not.ok)then
print *,msg
error stop
end if
end
real(dp) function determinant4(a) result(v)
real(dp),intent(in)::a(4,4)
real(dp)::q(4,4)
integer::ip(4),info,i
interface
subroutine dgetrf(m,n,a,lda,ipiv,info)
integer m,n,lda,ipiv(*),info
double precision a(lda,*)
end subroutine
end interface
q=a
call dgetrf(4,4,q,4,ip,info)
v=1d0
do i=1,4
v=v*q(i,i)
if(ip(i)/=i)v=-v
end do
end
end program
