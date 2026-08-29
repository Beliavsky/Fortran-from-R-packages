! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_stats
   use r_compat, only: dp, pchisq
   implicit none
   private
   public :: lrt, cov_to_cor
contains
   function lrt(loglik1,df1,loglik2,df2) result(out)
      real(dp),intent(in)::loglik1,loglik2
      integer,intent(in)::df1,df2
      real(dp)::out(2)
      integer::d
      out(1)=2.0_dp*abs(loglik2-loglik1)
      d=abs(df2-df1)
      if(d==0) then
         out(2)=merge(1.0_dp,0.0_dp,out(1)==0.0_dp)
      else
         out(2)=1.0_dp-pchisq(out(1),real(d,dp))
      end if
   end function lrt

   function cov_to_cor(cov) result(cor)
      real(dp),intent(in)::cov(:,:)
      real(dp)::cor(size(cov,1),size(cov,2)),den
      integer::i,j
      if(size(cov,1)/=size(cov,2))error stop 'cov_to_cor: matrix must be square'
      do i=1,size(cov,1)
         do j=1,size(cov,2)
            den=sqrt(max(0.0_dp,cov(i,i)*cov(j,j)))
            if(den>0.0_dp)then
            cor(i,j)=cov(i,j)/den
            else
            cor(i,j)=0.0_dp
            end if
         end do
      end do
      do i=1,size(cov,1)
      cor(i,i)=1.0_dp
      end do
   end function cov_to_cor
end module matrixdist_stats
