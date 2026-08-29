! SPDX-License-Identifier: GPL-3.0-only
module evd_bvpot
   use r_compat, only : dp
   use evd_bivariate, only : abvlog, abvalog, abvhr, abvneglog, abvaneglog, &
      abvbilog, abvnegbilog, abvct, abvamix, hbvlog, hbvhr, hbvneglog, hbvbilog, hbvnegbilog, hbvct
   implicit none
   private
   public :: bvpot_censored_nll, bvpot_poisson_nll
contains

function bvpot_censored_nll(data, threshold, model, dep, scale, shape) result(nll)
   real(dp), intent(in) :: data(:,:), threshold(2), dep(:), scale(2), shape(2)
   character(len=*), intent(in) :: model
   real(dp) :: nll
   integer :: n,i,nn,thid
   real(dp) :: lambda(2), e(2), xx(2), jac(2), xt(2), v,v1,v2,v12,ld
   logical :: hi1,hi2
   n=size(data,1)
   nll=0.0_dp
   if(size(data,2)/=2 .or. any(scale<=0.0_dp)) then
   nll=1.0e100_dp
   return
   end if
   lambda(1)=real(count(data(:,1)>threshold(1)),dp)/real(n+1,dp)
   lambda(2)=real(count(data(:,2)>threshold(2)),dp)/real(n+1,dp)
   if(any(lambda<=0.0_dp) .or. any(lambda>=1.0_dp)) then
   nll=1.0e100_dp
   return
   end if
   nn=0
   do i=1,n
      hi1=data(i,1)>threshold(1)
      hi2=data(i,2)>threshold(2)
      if(.not.(hi1.or.hi2)) cycle
      nn=nn+1
      thid=merge(1,0,hi1)+2*merge(1,0,hi2)
      e=[max(0.0_dp,data(i,1)-threshold(1)),max(0.0_dp,data(i,2)-threshold(2))]
      call pot_transform(e(1),lambda(1),scale(1),shape(1),xx(1),jac(1),hi1)
      if(.not.hi1) then
      nll=1.0e100_dp
      return
      end if
      call pot_transform(e(2),lambda(2),scale(2),shape(2),xx(2),jac(2),hi2)
      if(.not.hi2) then
      nll=1.0e100_dp
      return
      end if
      call exponent_derivs(xx(1),xx(2),model,dep,v,v1,v2,v12)
      select case(thid)
      case(1)
         if(-v1<=0.0_dp .or. jac(1)<=0.0_dp) then
         nll=1.0e100_dp
         return
         end if
         ld=log(-v1)+log(jac(1))-v
      case(2)
         if(-v2<=0.0_dp .or. jac(2)<=0.0_dp) then
         nll=1.0e100_dp
         return
         end if
         ld=log(-v2)+log(jac(2))-v
      case(3)
         if(v1*v2-v12<=0.0_dp .or. any(jac<=0.0_dp)) then
         nll=1.0e100_dp
         return
         end if
         ld=log(v1*v2-v12)+log(jac(1))+log(jac(2))-v
      end select
      nll=nll-ld
   end do
   xt=-1.0_dp/log(1.0_dp-lambda)
   v=exponent_value(1.0_dp/xt(1),1.0_dp/xt(2),model,dep)
   nll=nll+real(n-nn,dp)*v
end function bvpot_censored_nll

function bvpot_poisson_nll(data, threshold, model, dep, scale, shape) result(nll)
   real(dp), intent(in) :: data(:,:), threshold(2), dep(:), scale(2), shape(2)
   character(len=*), intent(in) :: model
   real(dp) :: nll
   integer :: n,i,thid
   real(dp) :: lambda(2), rprob(2), e(2), xx(2), jacm(2), rad,w,h,utt(2),a2,ld
   real(dp), allocatable :: rank1(:),rank2(:)
   logical :: hi1,hi2
   n=size(data,1)
   nll=0.0_dp
   if(size(data,2)/=2 .or. any(scale<=0.0_dp)) then
   nll=1.0e100_dp
   return
   end if
   if(trim(model)=='alog' .or. trim(model)=='aneglog' .or. trim(model)=='amix') then
      nll=1.0e100_dp
      return
   end if
   lambda(1)=real(count(data(:,1)>threshold(1)),dp)/real(n+1,dp)
   lambda(2)=real(count(data(:,2)>threshold(2)),dp)/real(n+1,dp)
   if(any(lambda<=0.0_dp) .or. any(lambda>=1.0_dp)) then
   nll=1.0e100_dp
   return
   end if
   allocate(rank1(n),rank2(n))
   call upper_rank_probs(data(:,1),rank1)
   call upper_rank_probs(data(:,2),rank2)
   do i=1,n
      hi1=data(i,1)>threshold(1)
      hi2=data(i,2)>threshold(2)
      if(.not.(hi1.or.hi2)) cycle
      thid=merge(1,0,hi1)+2*merge(1,0,hi2)
      rprob=[rank1(i),rank2(i)]
      if(hi1) rprob(1)=lambda(1)
      if(hi2) rprob(2)=lambda(2)
      e=[max(0.0_dp,data(i,1)-threshold(1)),max(0.0_dp,data(i,2)-threshold(2))]
      call pot_pp_transform(e(1),rprob(1),lambda(1),scale(1),shape(1),xx(1),jacm(1),hi1)
      if(.not.hi1) then
      nll=1.0e100_dp
      return
      end if
      call pot_pp_transform(e(2),rprob(2),lambda(2),scale(2),shape(2),xx(2),jacm(2),hi2)
      if(.not.hi2) then
      nll=1.0e100_dp
      return
      end if
      rad=xx(1)+xx(2)
      if(rad<=0.0_dp) then
      nll=1.0e100_dp
      return
      end if
      w=xx(1)/rad
      h=max(model_h(w,model,dep),tiny(1.0_dp))
      select case(thid)
      case(1); ld=jacm(1)+log(h)-3.0_dp*log(rad)
      case(2); ld=jacm(2)+log(h)-3.0_dp*log(rad)
      case(3); ld=jacm(1)+jacm(2)+log(h)-3.0_dp*log(rad)
      end select
      nll=nll-ld
   end do
   utt=-1.0_dp/log(1.0_dp-lambda)
   nll=nll+exponent_value(1.0_dp/utt(1),1.0_dp/utt(2),model,dep)
end function bvpot_poisson_nll

pure subroutine pot_transform(excess,lambda,scale,shape,x,jac,ok)
   real(dp),intent(in)::excess,lambda,scale,shape
   real(dp),intent(out)::x,jac
   logical,intent(out)::ok
   real(dp)::z,t
   ok=.false.
   if(scale<=0.0_dp) return
   z=excess/scale
   if(abs(shape)<=epsilon(1.0_dp)**0.3_dp) then
   t=exp(-z)
   else
      if(1.0_dp+shape*z<=0.0_dp) return
      t=(1.0_dp+shape*z)**(-1.0_dp/shape)
   end if
   if(lambda*t>=1.0_dp) return
   x=-1.0_dp/log(1.0_dp-lambda*t)
   jac=lambda*x*x*t**(1.0_dp+shape)/(scale*(1.0_dp-lambda*t))
   ok=jac>0.0_dp
end subroutine pot_transform

pure subroutine pot_pp_transform(excess,rprob,p,scale,shape,x,logjac,ok)
   real(dp),intent(in)::excess,rprob,p,scale,shape
   real(dp),intent(out)::x,logjac
   logical,intent(out)::ok
   real(dp)::z,t,om
   ok=.false.
   if(scale<=0.0_dp .or. rprob<=0.0_dp .or. rprob>=1.0_dp .or. p<=0.0_dp) return
   z=excess/scale
   if(abs(shape)<=epsilon(1.0_dp)**0.3_dp) then
   t=exp(-z)
   else
      if(1.0_dp+shape*z<=0.0_dp) return
      t=(1.0_dp+shape*z)**(-1.0_dp/shape)
   end if
   if(rprob*t>=1.0_dp) return
   x=-1.0_dp/log(1.0_dp-rprob*t)
   om=max(1.0_dp-exp(-1.0_dp/x),tiny(1.0_dp))
   logjac=2.0_dp*log(x)+1.0_dp/x+(1.0_dp+shape)*log(om)-log(scale)-shape*log(p)
   ok=.true.
end subroutine pot_pp_transform

pure function exponent_value(x1,x2,model,dep) result(v)
   real(dp),intent(in)::x1,x2,dep(:)
   character(len=*),intent(in)::model
   real(dp)::v,s,w,a
   s=x1+x2
   if(s<=0.0_dp) then
   v=0.0_dp
   return
   end if
   w=x1/s
   a=model_a(w,model,dep)
   v=s*a
end function exponent_value

pure subroutine exponent_derivs(x1,x2,model,dep,v,v1,v2,v12)
   real(dp),intent(in)::x1,x2,dep(:)
   character(len=*),intent(in)::model
   real(dp),intent(out)::v,v1,v2,v12
   real(dp)::s,w,a,a1,a2
   real(dp) :: z1,z2,vx1,vx2,vx12
   ! x1,x2 here are unit-Frechet coordinates z.  The Pickands
   ! representation is expressed in reciprocal coordinates 1/z.
   z1=x1
   z2=x2
   s=1.0_dp/z1+1.0_dp/z2
   w=(1.0_dp/z1)/s
   call a_derivs(w,model,dep,a,a1,a2)
   v=s*a
   vx1=a+(1.0_dp-w)*a1
   vx2=a-w*a1
   vx12=-w*(1.0_dp-w)*a2/s
   v1=-vx1/(z1*z1)
   v2=-vx2/(z2*z2)
   v12=vx12/(z1*z1*z2*z2)
end subroutine exponent_derivs

pure subroutine a_derivs(w,model,dep,a,a1,a2)
   real(dp),intent(in)::w,dep(:)
   character(len=*),intent(in)::model
   real(dp),intent(out)::a,a1,a2
   real(dp)::h,wm,wp,am,ap
   h=max(1.0e-5_dp,epsilon(1.0_dp)**0.25_dp)
   h=min(h,0.25_dp*max(min(w,1.0_dp-w),h))
   wm=max(1.0e-9_dp,w-h)
   wp=min(1.0_dp-1.0e-9_dp,w+h)
   a=model_a(w,model,dep)
   am=model_a(wm,model,dep)
   ap=model_a(wp,model,dep)
   a1=(ap-am)/(wp-wm)
   a2=2.0_dp*((ap-a)/(wp-w)-(a-am)/(w-wm))/(wp-wm)
end subroutine a_derivs
pure function a_second(w,model,dep) result(a2)
   real(dp),intent(in)::w,dep(:)
   character(len=*),intent(in)::model
   real(dp)::a2,a,a1
   call a_derivs(min(max(w,1.0e-8_dp),1.0_dp-1.0e-8_dp),model,dep,a,a1,a2)
end function a_second

pure function model_a(w,model,p) result(a)
   real(dp),intent(in)::w,p(:)
   character(len=*),intent(in)::model
   real(dp)::a
   select case(trim(model))
   case('log'); a=abvlog(w,p(1))
   case('alog'); a=abvalog(w,p(1),p(2),p(3))
   case('hr'); a=abvhr(w,p(1))
   case('neglog'); a=abvneglog(w,p(1))
   case('aneglog'); a=abvaneglog(w,p(1),p(2),p(3))
   case('bilog'); a=abvbilog(w,p(1),p(2))
   case('negbilog'); a=abvnegbilog(w,p(1),p(2))
   case('ct'); a=abvct(w,p(1),p(2))
   case('amix'); a=abvamix(w,p(1),p(2))
   case default; a=huge(1.0_dp)
   end select
end function model_a


pure function model_h(w,model,p) result(h)
   real(dp),intent(in)::w,p(:)
   character(len=*),intent(in)::model
   real(dp)::h
   select case(trim(model))
   case('log'); h=hbvlog(w,p(1))
   case('hr'); h=hbvhr(w,p(1))
   case('neglog'); h=hbvneglog(w,p(1))
   case('bilog'); h=hbvbilog(w,p(1),p(2))
   case('negbilog'); h=hbvnegbilog(w,p(1),p(2))
   case('ct'); h=hbvct(w,p(1),p(2))
   case default; h=0.0_dp
   end select
end function model_h

subroutine upper_rank_probs(x,p)
   real(dp),intent(in)::x(:)
   real(dp),intent(out)::p(size(x))
   integer::i,j
   real(dp)::rank
   do i=1,size(x)
      rank=1.0_dp
      do j=1,size(x)
         if(j==i) cycle
         if(x(j)<x(i)) rank=rank+1.0_dp
         if(x(j)==x(i)) rank=rank+0.5_dp
      end do
      p(i)=1.0_dp-rank/real(size(x)+1,dp)
      p(i)=min(max(p(i),tiny(1.0_dp)),1.0_dp-epsilon(1.0_dp))
   end do
end subroutine upper_rank_probs
end module evd_bvpot
