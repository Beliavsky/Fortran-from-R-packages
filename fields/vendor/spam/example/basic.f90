program basic
use spam
implicit none
real(dp)::x(3,3),b(3);real(dp),allocatable::sol(:)
type(csr_matrix)::a;type(spam_chol)::f
x=reshape([4d0,1d0,0d0,1d0,3d0,1d0,0d0,1d0,2d0],[3,3]);a=csr_from_dense(x)
f=spam_chol_factor(a);b=[1d0,2d0,3d0];sol=spam_solve(f,b)
print '(a,3f12.6)','solution: ',sol
end program
