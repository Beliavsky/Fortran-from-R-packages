! SPDX-License-Identifier: GPL-3.0-only
program test_sparse_ldlt
   use scs_kinds, only : dp, i4
   use scs_types, only : scs_matrix
   use scs_sparse, only : dense_to_csc, dense_upper_to_csc
   use scs_ldlt, only : scs_ldlt_factor
   implicit none
   type(scs_matrix) :: a, p
   type(scs_ldlt_factor) :: fac
   real(dp) :: ad(5,4), pd(4,4), diag(9), rhs(9), rhs0(9), res(9)
   logical :: ok

   ad = 0.0_dp
   ad(1,1)=1.0_dp; ad(1,2)=0.5_dp
   ad(2,2)=-2.0_dp; ad(2,3)=1.0_dp
   ad(3,1)=0.25_dp; ad(3,4)=1.5_dp
   ad(4,3)=-0.75_dp
   ad(5,2)=1.25_dp; ad(5,4)=-0.5_dp
   pd = 0.0_dp
   pd(1,1)=2.0_dp; pd(2,2)=1.0_dp; pd(3,3)=3.0_dp; pd(4,4)=0.5_dp
   pd(1,2)=0.2_dp; pd(2,3)=-0.1_dp; pd(3,4)=0.15_dp
   call dense_to_csc(ad,a)
   call dense_upper_to_csc(pd,p)

   diag(1:4)=1.0e-3_dp
   diag(5:9)=0.2_dp
   call fac%factorize(a,p,.true.,diag,ok)
   if(.not.ok) error stop 'sparse LDL initial factorization failed'
   if(fac%symbolic_analyses/=1_i4 .or. fac%factorizations/=1_i4) error stop 'bad initial factorization counters'

   rhs0=[1.0_dp,-2.0_dp,0.5_dp,3.0_dp,-1.0_dp,2.0_dp,0.25_dp,-0.5_dp,1.5_dp]
   rhs=rhs0
   call fac%solve(rhs)
   call sym_upper_mv(fac%kkt,rhs,res)
   if(maxval(abs(res-rhs0))>5.0e-11_dp) then
      print *, 'residual = ', maxval(abs(res-rhs0))
      error stop 'sparse LDL solve residual'
   end if

   diag(1:4)=2.0e-3_dp
   diag(5:9)=0.35_dp
   call fac%factorize(a,p,.true.,diag,ok)
   if(.not.ok) error stop 'sparse LDL refactorization failed'
   if(fac%symbolic_analyses/=1_i4 .or. fac%factorizations/=2_i4) error stop 'symbolic analysis was not reused'
   if(fac%kkt_nnz/=int(size(fac%kkt%x),i4)) error stop 'bad KKT nnz statistic'
   if(fac%factor_nnz/=int(size(fac%lx),i4)) error stop 'bad factor nnz statistic'

   rhs=rhs0
   call fac%solve(rhs)
   call sym_upper_mv(fac%kkt,rhs,res)
   if(maxval(abs(res-rhs0))>5.0e-11_dp) error stop 'sparse LDL refactorized solve residual'

   print '(a)', 'Sparse QDLDL backend tests passed.'
contains
   subroutine sym_upper_mv(s,x,y)
      type(scs_matrix),intent(in)::s
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      integer(i4)::c,k,r
      y=0.0_dp
      do c=1_i4,s%n
         do k=s%p(c),s%p(c+1_i4)-1_i4
            r=s%i(k)
            y(r)=y(r)+s%x(k)*x(c)
            if(r/=c)y(c)=y(c)+s%x(k)*x(r)
         end do
      end do
   end subroutine sym_upper_mv
end program test_sparse_ldlt
