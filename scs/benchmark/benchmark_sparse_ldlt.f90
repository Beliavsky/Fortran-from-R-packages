! SPDX-License-Identifier: GPL-3.0-only
program benchmark_sparse_ldlt
   use scs_kinds, only : dp
   use scs_types, only : scs_matrix
   use scs_sparse, only : dense_to_csc, dense_upper_to_csc
   use scs_ldlt, only : scs_ldlt_factor
   implicit none
   integer, parameter :: nv=300, mc=300
   type(scs_matrix) :: a,p
   type(scs_ldlt_factor) :: fac
   real(dp), allocatable :: ad(:,:),pd(:,:),diag(:),rhs(:),rhs0(:),k(:,:),l(:,:),d(:)
   real(dp) :: t0,t1,tsparse,tdense,sparse_mb,dense_mb
   integer :: i,j,r,kk,nm
   logical :: ok

   nm=nv+mc
   allocate(ad(mc,nv),pd(nv,nv),diag(nm),rhs(nm),rhs0(nm))
   ad=0.0_dp;pd=0.0_dp
   do j=1,nv
      pd(j,j)=1.0_dp+0.01_dp*real(mod(j,7),dp)
      do r=max(1,j-2),min(mc,j+2)
         ad(r,j)=0.25_dp+0.01_dp*real(mod(3*r+5*j,17),dp)
      end do
   end do
   call dense_to_csc(ad,a)
   call dense_upper_to_csc(pd,p)
   diag(1:nv)=1.0e-3_dp;diag(nv+1:nm)=0.2_dp

   call cpu_time(t0)
   call fac%factorize(a,p,.true.,diag,ok)
   call cpu_time(t1)
   if(.not.ok)error stop 'sparse factorization failed'
   tsparse=t1-t0
   rhs0=[(sin(0.1_dp*real(i,dp)),i=1,nm)];rhs=rhs0
   call fac%solve(rhs)

   allocate(k(nm,nm),source=0.0_dp)
   do j=1,fac%kkt%n
      do kk=fac%kkt%p(j),fac%kkt%p(j+1)-1
         r=fac%kkt%i(kk);k(r,j)=k(r,j)+fac%kkt%x(kk)
         if(r/=j)k(j,r)=k(j,r)+fac%kkt%x(kk)
      end do
   end do
   allocate(l(nm,nm),d(nm));l=0.0_dp
   do i=1,nm;l(i,i)=1.0_dp;end do
   call cpu_time(t0)
   call dense_ldlt(k,l,d,ok)
   call cpu_time(t1)
   if(.not.ok)error stop 'dense factorization failed'
   tdense=t1-t0

   sparse_mb=real(size(fac%lx),dp)*12.0_dp/1048576.0_dp + real(nm,dp)*8.0_dp/1048576.0_dp
   dense_mb=real(nm,dp)*real(nm,dp)*8.0_dp/1048576.0_dp
   print '(a,i0)', 'KKT dimension: ',nm
   print '(a,i0)', 'KKT nnz (upper): ',fac%kkt_nnz
   print '(a,i0)', 'L nnz (strict lower): ',fac%factor_nnz
   print '(a,f10.6,a)', 'Sparse QDLDL factor time: ',tsparse,' s'
   print '(a,f10.6,a)', 'Dense LDL factor time:     ',tdense,' s'
   if(tsparse>0.0_dp)print '(a,f10.2,a)','Dense/sparse time ratio:  ',tdense/tsparse,'x'
   print '(a,f10.3,a)', 'Approx sparse factor payload: ',sparse_mb,' MiB'
   print '(a,f10.3,a)', 'Dense L payload:             ',dense_mb,' MiB'
contains
   subroutine dense_ldlt(a,l,d,good)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(inout)::l(:,:)
      real(dp),intent(out)::d(:)
      logical,intent(out)::good
      integer::jj,rr,k0
      real(dp)::s
      good=.true.
      do jj=1,size(a,1)
         s=a(jj,jj)
         do k0=1,jj-1;s=s-l(jj,k0)*l(jj,k0)*d(k0);end do
         d(jj)=s
         if(abs(s)<=tiny(1.0_dp))then;good=.false.;return;end if
         do rr=jj+1,size(a,1)
            s=a(rr,jj)
            do k0=1,jj-1;s=s-l(rr,k0)*l(jj,k0)*d(k0);end do
            l(rr,jj)=s/d(jj)
         end do
      end do
   end subroutine dense_ldlt
end program benchmark_sparse_ldlt
