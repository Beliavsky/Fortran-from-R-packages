module lavaan_categorical_analytic
   use lavaan_kinds, only : dp
   use lavaan_categorical, only : categorical_stats_result
   use lavaan_ordinal, only : ordinal_thresholds, polychoric_matrix, bvn_rectangle
   use lavaan_linalg, only : inverse_general
   implicit none
   private
   public :: categorical_wls_statistics_analytic
contains
   subroutine categorical_wls_statistics_analytic(data,result)
      integer,intent(in)::data(:,:)
      type(categorical_stats_result),intent(out)::result
      real(dp),allocatable::theta(:),g(:,:),gp(:,:),gm(:,:),gbar(:),a(:,:),ainv(:,:),b(:,:),d(:),ridge(:,:),winv(:,:)
      integer::n,p,q,r,j,info
      real(dp)::h
      n=size(data,1)
      p=size(data,2)
      if(n<5 .or. p<2 .or. minval(data)<1) then
      result%status=-1
      return
      end if
      call setup_stats(data,result,theta,info)
      if(info/=0) then
      result%status=info
      return
      end if
      q=size(theta)
      allocate(g(n,q))
      call estimating_scores(data,result,theta,g,info)
      if(info/=0) then
      result%status=info
      return
      end if
      gbar=sum(g,dim=1)/real(n,dp)
      allocate(a(q,q),gp(n,q),gm(n,q))
      a=0.0_dp
      do j=1,q
         h=2.0e-5_dp*max(1.0_dp,abs(theta(j)))
         block
            real(dp),allocatable::tp(:),tm(:)
            tp=theta
            tm=theta
            tp(j)=tp(j)+h
            tm(j)=tm(j)-h
            if(j>size(result%thresholds)) then
               tp(j)=min(0.998_dp,tp(j))
               tm(j)=max(-0.998_dp,tm(j))
               h=0.5_dp*(tp(j)-tm(j))
            end if
            call estimating_scores(data,result,tp,gp,info)
            if(info/=0) then
            result%status=info
            return
            end if
            call estimating_scores(data,result,tm,gm,info)
            if(info/=0) then
            result%status=info
            return
            end if
         end block
         a(:,j)=(sum(gp,dim=1)-sum(gm,dim=1))/(2.0_dp*h*real(n,dp))
      end do
      call inverse_general(a,ainv,info)
      if(info/=0) then
      result%status=200+info
      return
      end if
      allocate(b(q,q))
      b=0.0_dp
      do r=1,n
         d=g(r,:)-gbar
         b=b+spread(d,2,q)*spread(d,1,q)
      end do
      b=b/real(n,dp)
      result%gamma=matmul(ainv,matmul(b,transpose(ainv)))
      result%gamma=0.5_dp*(result%gamma+transpose(result%gamma))
      ridge=result%gamma
      do j=1,q
      ridge(j,j)=ridge(j,j)+1.0e-9_dp*max(1.0_dp,abs(ridge(j,j)))
      end do
      call inverse_general(ridge,winv,info)
      allocate(result%weight(q,q),result%dwls_weight(q))
      result%weight=0.0_dp
      result%dwls_weight=0.0_dp
      if(info==0) then
         result%weight=winv
      else
         do j=1,q
         if(result%gamma(j,j)>1.0e-12_dp) result%weight(j,j)=1.0_dp/result%gamma(j,j)
         end do
      end if
      do j=1,q
      result%dwls_weight(j)=result%weight(j,j)
      end do
      result%n_jackknife=0
      result%status=0
   end subroutine categorical_wls_statistics_analytic

   subroutine setup_stats(data,result,theta,info)
      integer,intent(in)::data(:,:)
      type(categorical_stats_result),intent(inout)::result
      real(dp),allocatable,intent(out)::theta(:)
      integer,intent(out)::info
      integer::p,j,i,pos,nth,ncor,k
      integer,allocatable::cnt(:)
      real(dp),allocatable::th(:)
      p=size(data,2)
      allocate(result%ncat(p),result%threshold_offset(p+1))
      result%threshold_offset(1)=1
      nth=0
      info=0
      do j=1,p
         result%ncat(j)=maxval(data(:,j))
         if(result%ncat(j)<2 .or. minval(data(:,j))<1) then
         info=j
         return
         end if
         nth=nth+result%ncat(j)-1
         result%threshold_offset(j+1)=nth+1
      end do
      allocate(result%thresholds(nth))
      pos=1
      do j=1,p
         allocate(cnt(result%ncat(j)))
         cnt=0
         do i=1,size(data,1)
         cnt(data(i,j))=cnt(data(i,j))+1
         end do
         if(any(cnt==0)) then
         info=100+j
         return
         end if
         th=ordinal_thresholds(cnt)
         result%thresholds(pos:pos+size(th)-1)=th
         pos=pos+size(th)
         deallocate(cnt)
      end do
      call polychoric_matrix(data,result%correlation,info)
      if(info/=0) return
      ncor=p*(p-1)/2
      allocate(theta(nth+ncor))
      theta(1:nth)=result%thresholds
      k=nth
      do j=1,p-1
      do i=j+1,p
      k=k+1
      theta(k)=result%correlation(i,j)
      end do
      end do
      result%stats=theta
   end subroutine setup_stats

   subroutine estimating_scores(data,result,theta,g,info)
      integer,intent(in)::data(:,:)
      type(categorical_stats_result),intent(in)::result
      real(dp),intent(in)::theta(:)
      real(dp),intent(out)::g(:,:)
      integer,intent(out)::info
      integer::n,p,nth,r,j,i,k,c1,c2,tidx
      real(dp)::rho,h,lp,lo1,hi1,lo2,hi2,prp,prm
      n=size(data,1)
      p=size(data,2)
      nth=size(result%thresholds)
      g=0.0_dp
      info=0
      ! Marginal threshold estimating equations: I(Y <= c) - Phi(tau_c).
      do j=1,p
         do c1=1,result%ncat(j)-1
            tidx=result%threshold_offset(j)+c1-1
            do r=1,n
               g(r,tidx)=merge(1.0_dp,0.0_dp,data(r,j)<=c1)-normal_cdf(theta(tidx))
            end do
         end do
      end do
      k=nth
      do j=1,p-1
         do i=j+1,p
            k=k+1
            rho=max(-0.998_dp,min(0.998_dp,theta(k)))
            h=2.0e-5_dp
            do r=1,n
               c1=data(r,j)
               c2=data(r,i)
               call category_bounds(j,c1,result,theta,lo1,hi1)
               call category_bounds(i,c2,result,theta,lo2,hi2)
               prp=bvn_rectangle(lo1,hi1,lo2,hi2,min(0.999_dp,rho+h))
               prm=bvn_rectangle(lo1,hi1,lo2,hi2,max(-0.999_dp,rho-h))
               lp=(log(max(prp,1.0e-300_dp))-log(max(prm,1.0e-300_dp)))/(2.0_dp*h)
               g(r,k)=lp
            end do
         end do
      end do
   end subroutine estimating_scores

   subroutine category_bounds(j,c,result,theta,lo,hi)
      integer,intent(in)::j,c
      type(categorical_stats_result),intent(in)::result
      real(dp),intent(in)::theta(:)
      real(dp),intent(out)::lo,hi
      integer::base
      base=result%threshold_offset(j)
      if(c==1) then
      lo=-huge(1.0_dp)
      else
      lo=theta(base+c-2)
      end if
      if(c==result%ncat(j)) then
      hi=huge(1.0_dp)
      else
      hi=theta(base+c-1)
      end if
   end subroutine category_bounds

   pure function normal_cdf(x) result(p)
      real(dp),intent(in)::x
      real(dp)::p
      if(x>8.0_dp) then
      p=1.0_dp
      else if(x< -8.0_dp) then
      p=0.0_dp
      else
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
      end if
   end function normal_cdf
end module lavaan_categorical_analytic
