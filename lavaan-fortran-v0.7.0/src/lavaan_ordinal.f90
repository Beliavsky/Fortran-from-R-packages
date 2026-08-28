module lavaan_ordinal
   use lavaan_kinds, only : dp
   use pbivnorm_mod, only : pbivnorm
   implicit none
   private
   public :: ordinal_thresholds, polychoric_table, polychoric_matrix, bvn_rectangle, normal_quantile
contains
   pure function normal_cdf(x) result(p)
      real(dp),intent(in)::x
      real(dp)::p
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   function normal_quantile(p) result(x)
      real(dp),intent(in)::p
      real(dp)::x
      real(dp),parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp
      real(dp),parameter :: a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp
      real(dp),parameter :: a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
      real(dp),parameter :: b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp
      real(dp),parameter :: b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp
      real(dp),parameter :: c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp
      real(dp),parameter :: c3=-2.400758277161838_dp,c4=-2.549732539343734_dp
      real(dp),parameter :: c5=4.374664141464968_dp,c6=2.938163982698783_dp
      real(dp),parameter :: d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp
      real(dp),parameter :: d3=2.445134137142996_dp,d4=3.754408661907416_dp
      real(dp)::q,r,phi
      if(p<=0) then
      x=-huge(1.0_dp)
      return
      else if(p>=1) then
      x=huge(1.0_dp)
      return
      end if
      if(p<0.02425_dp) then
         q=sqrt(-2*log(p))
         x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      else if(p>0.97575_dp) then
         q=sqrt(-2*log(1-p))
         x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
      else
         q=p-0.5_dp
         r=q*q
         x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      end if
      phi=exp(-0.5_dp*x*x)/sqrt(2*acos(-1.0_dp))
      x=x-(normal_cdf(x)-p)/phi
   end function normal_quantile

   function ordinal_thresholds(counts) result(th)
      integer,intent(in)::counts(:)
      real(dp),allocatable::th(:)
      integer::k,n,cum
      n=sum(counts)
      allocate(th(max(0,size(counts)-1)))
      cum=0
      do k=1,size(th)
      cum=cum+counts(k)
      th(k)=normal_quantile(real(cum,dp)/real(n,dp))
      end do
   end function ordinal_thresholds

   function bvn_cdf(x,y,rho) result(p)
      real(dp),intent(in)::x,y,rho
      real(dp)::p
      if(x < -8.0_dp .or. y < -8.0_dp) then
         p=0.0_dp
      else if(x > 8.0_dp .and. y > 8.0_dp) then
         p=1.0_dp
      else if(x > 8.0_dp) then
         p=normal_cdf(y)
      else if(y > 8.0_dp) then
         p=normal_cdf(x)
      else
         p=pbivnorm(x,y,rho)
      end if
   end function bvn_cdf

   function bvn_rectangle(l1,u1,l2,u2,rho) result(p)
      real(dp),intent(in)::l1,u1,l2,u2,rho
      real(dp)::p
      p=bvn_cdf(u1,u2,rho)-bvn_cdf(l1,u2,rho)-bvn_cdf(u1,l2,rho)+bvn_cdf(l1,l2,rho)
      p=max(p,1.0e-300_dp)
   end function bvn_rectangle

   subroutine polychoric_table(tab,rho,threshold1,threshold2,loglik)
      integer,intent(in)::tab(:,:)
      real(dp),intent(out)::rho,loglik
      real(dp),allocatable,intent(out)::threshold1(:),threshold2(:)
      integer,allocatable::r(:),c(:)
      real(dp)::a,b,x1,x2,f1,f2,gr
      integer::iter
      allocate(r(size(tab,1)),c(size(tab,2)))
      r=sum(tab,dim=2)
      c=sum(tab,dim=1)
      threshold1=ordinal_thresholds(r)
      threshold2=ordinal_thresholds(c)
      gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      a=-0.999_dp
      b=0.999_dp
      x1=b-gr*(b-a)
      x2=a+gr*(b-a)
      f1=negll(x1)
      f2=negll(x2)
      do iter=1,100
         if(abs(b-a)<1.0e-8_dp) exit
         if(f1>f2) then
         a=x1
         x1=x2
         f1=f2
         x2=a+gr*(b-a)
         f2=negll(x2)
         else
         b=x2
         x2=x1
         f2=f1
         x1=b-gr*(b-a)
         f1=negll(x1)
         end if
      end do
      rho=0.5_dp*(a+b)
      loglik=-negll(rho)
   contains
      function negll(rr) result(v)
         real(dp),intent(in)::rr
         real(dp)::v,l1,u1,l2,u2,pr
         integer::ii,jj
         v=0
         do ii=1,size(tab,1)
            if(ii==1) then
            l1=-huge(1.0_dp)
            else
            l1=threshold1(ii-1)
            end if
            if(ii==size(tab,1)) then
            u1=huge(1.0_dp)
            else
            u1=threshold1(ii)
            end if
            do jj=1,size(tab,2)
               if(jj==1) then
               l2=-huge(1.0_dp)
               else
               l2=threshold2(jj-1)
               end if
               if(jj==size(tab,2)) then
               u2=huge(1.0_dp)
               else
               u2=threshold2(jj)
               end if
               pr=bvn_rectangle(l1,u1,l2,u2,rr)
               if(tab(ii,jj)>0) v=v-real(tab(ii,jj),dp)*log(pr)
            end do
         end do
      end function negll
   end subroutine polychoric_table

   subroutine polychoric_matrix(data, cor, info)
      integer, intent(in) :: data(:, :)
      real(dp), allocatable, intent(out) :: cor(:, :)
      integer, intent(out) :: info
      integer :: p, i, j, a, ni, nj
      integer, allocatable :: tab(:, :)
      real(dp), allocatable :: t1(:), t2(:)
      real(dp) :: rho, ll
      p=size(data,2)
      allocate(cor(p,p))
      cor=0.0_dp
      info=0
      do i=1,p
      cor(i,i)=1.0_dp
      end do
      do j=1,p-1
         nj=maxval(data(:,j))
         if(minval(data(:,j))<1 .or. nj<2) then
         info=j
         return
         end if
         do i=j+1,p
            ni=maxval(data(:,i))
            if(minval(data(:,i))<1 .or. ni<2) then
            info=i
            return
            end if
            allocate(tab(nj,ni))
            tab=0
            do a=1,size(data,1)
            tab(data(a,j),data(a,i))=tab(data(a,j),data(a,i))+1
            end do
            call polychoric_table(tab,rho,t1,t2,ll)
            cor(j,i)=rho
            cor(i,j)=rho
            deallocate(tab,t1,t2)
         end do
      end do
   end subroutine polychoric_matrix
end module lavaan_ordinal
