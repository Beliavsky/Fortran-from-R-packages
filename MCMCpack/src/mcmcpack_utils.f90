! SPDX-License-Identifier: GPL-3.0-only
module mcmcpack_utils
   use mcmcpack_kinds, only : dp
   use mcmcpack_math, only : logsumexp
   implicit none
   private
   type, public :: procrustes_result
      real(dp), allocatable :: x_new(:,:), rotation(:,:), translation(:)
      real(dp) :: scale = 1.0_dp
      integer :: status = 0
   end type procrustes_result
   type, public :: waic_result
      real(dp) :: waic=0.0_dp,waic2=0.0_dp,elpd_waic=0.0_dp,p_waic=0.0_dp
      real(dp) :: p_waic1=0.0_dp,elpd_loo=0.0_dp,p_loo=0.0_dp
      real(dp) :: waic_ci(2)=0.0_dp
      real(dp), allocatable :: pointwise(:,:), total(:), se(:)
   end type waic_result
   public :: vech,xpnd,procrustes,waic
contains
   function vech(x) result(v)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable::v(:)
      integer::n,i,j,k
      if(size(x,1)/=size(x,2))then;allocate(v(0));return;end if
      n=size(x,1);allocate(v(n*(n+1)/2));k=0
      do j=1,n;do i=j,n;k=k+1;v(k)=x(i,j);end do;end do
   end function vech

   function xpnd(v,nrow) result(x)
      real(dp),intent(in)::v(:)
      integer,intent(in),optional::nrow
      real(dp),allocatable::x(:,:)
      integer::n,i,j,k
      if(present(nrow))then;n=nrow;else;n=nint((-1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(v),dp)))/2.0_dp);end if
      allocate(x(n,n));x=0.0_dp;k=0
      do j=1,n;do i=j,n;k=k+1;x(i,j)=v(mod(k-1,size(v))+1);x(j,i)=x(i,j);end do;end do
   end function xpnd

   subroutine jacobi_sym(a,eval,evec,info)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::eval(size(a,1)),evec(size(a,1),size(a,1))
      integer,intent(out)::info
      real(dp)::b(size(a,1),size(a,1)),app,aqq,apq,tau,t,c,s,bip,biq,vip,viq
      integer::n,p,q,i,iter,maxiter
      n=size(a,1);b=a;evec=0.0_dp;do i=1,n;evec(i,i)=1.0_dp;end do
      maxiter=100*n*n;info=0
      do iter=1,maxiter
         p=1;q=min(2,n);apq=0.0_dp
         do i=1,n
            if(i<n)then
               if(maxval(abs(b(i,i+1:n)))>abs(apq))then;q=i+maxloc(abs(b(i,i+1:n)),dim=1);p=i;apq=b(p,q);end if
            end if
         end do
         if(abs(apq)<100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(b))))exit
         app=b(p,p);aqq=b(q,q);tau=(aqq-app)/(2.0_dp*apq)
         if(tau>=0.0_dp)then;t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau));else;t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau));end if
         c=1.0_dp/sqrt(1.0_dp+t*t);s=t*c
         do i=1,n
            if(i/=p.and.i/=q)then
               bip=b(i,p);biq=b(i,q);b(i,p)=c*bip-s*biq;b(p,i)=b(i,p);b(i,q)=s*bip+c*biq;b(q,i)=b(i,q)
            end if
            vip=evec(i,p);viq=evec(i,q);evec(i,p)=c*vip-s*viq;evec(i,q)=s*vip+c*viq
         end do
         b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq;b(p,q)=0.0_dp;b(q,p)=0.0_dp
      end do
      if(iter>maxiter)info=1
      do i=1,n;eval(i)=b(i,i);end do
      call sort_eigen_desc(eval,evec)
   end subroutine jacobi_sym

   subroutine sort_eigen_desc(eval,evec)
      real(dp),intent(inout)::eval(:),evec(:,:)
      integer::i,j,k,n
      real(dp)::tmp,col(size(eval))
      n=size(eval)
      do i=1,n-1
         k=i;do j=i+1,n;if(eval(j)>eval(k))k=j;end do
         if(k/=i)then;tmp=eval(i);eval(i)=eval(k);eval(k)=tmp;col=evec(:,i);evec(:,i)=evec(:,k);evec(:,k)=col;end if
      end do
   end subroutine sort_eigen_desc

   function procrustes(x,xstar,translation,dilation) result(res)
      real(dp),intent(in)::x(:,:),xstar(:,:)
      logical,intent(in),optional::translation,dilation
      type(procrustes_result)::res
      logical::tr,di
      integer::n,m,j,info
      real(dp),allocatable::xc(:,:),yc(:,:),c(:,:),ata(:,:),v(:,:),u(:,:),eval(:),sing(:),r(:,:)
      real(dp)::mx(size(x,2)),my(size(x,2)),num,den
      tr=.false.;di=.false.;if(present(translation))tr=translation;if(present(dilation))di=dilation
      if(any(shape(x)/=shape(xstar)))then;res%status=1;return;end if
      n=size(x,1);m=size(x,2);allocate(xc(n,m),yc(n,m),c(m,m),ata(m,m),v(m,m),u(m,m),eval(m),sing(m),r(m,m))
      mx=0.0_dp;my=0.0_dp;if(tr)then;mx=sum(x,dim=1)/real(n,dp);my=sum(xstar,dim=1)/real(n,dp);end if
      do j=1,m;xc(:,j)=x(:,j)-mx(j);yc(:,j)=xstar(:,j)-my(j);end do
      c=matmul(transpose(yc),xc);ata=matmul(transpose(c),c);call jacobi_sym(ata,eval,v,info)
      if(info/=0)then;res%status=2;return;end if
      sing=sqrt(max(eval,0.0_dp));u=0.0_dp
      do j=1,m
         if(sing(j)>sqrt(epsilon(1.0_dp)))u(:,j)=matmul(c,v(:,j))/sing(j)
      end do
      ! Complete only degenerate columns with coordinate basis, then normalize.
      do j=1,m
         if(sum(u(:,j)*u(:,j))<=epsilon(1.0_dp))u(j,j)=1.0_dp
         u(:,j)=u(:,j)/sqrt(max(dot_product(u(:,j),u(:,j)),tiny(1.0_dp)))
      end do
      r=matmul(v,transpose(u));res%scale=1.0_dp
      if(di)then
         num=sum(matmul(xc,r)*yc);den=sum(xc*xc);if(den>0.0_dp)res%scale=num/den
      end if
      allocate(res%translation(m));res%translation=0.0_dp
      if(tr)res%translation=my-res%scale*matmul(mx,r)
      allocate(res%rotation(m,m),res%x_new(n,m));res%rotation=r
      res%x_new=res%scale*matmul(x,r)
      do j=1,m;res%x_new(:,j)=res%x_new(:,j)+res%translation(j);end do
   end function procrustes

   function waic(log_lik) result(res)
      real(dp),intent(in)::log_lik(:,:)
      type(waic_result)::res
      integer::s,n,i,j
      real(dp)::lpd,pw,pw1,elw,eloo,ploo,mx,m1,m2,var,wmean,weighted,wsum
      real(dp),allocatable::w(:),vals(:)
      s=size(log_lik,1);n=size(log_lik,2)
      allocate(res%pointwise(n,8),res%total(8),res%se(8),w(s),vals(n));res%pointwise=0.0_dp
      do j=1,n
         mx=maxval(log_lik(:,j));lpd=mx+log(sum(exp(log_lik(:,j)-mx))/real(s,dp))
         m1=sum(log_lik(:,j))/real(s,dp);pw=sum((log_lik(:,j)-m1)**2)/real(max(1,s-1),dp)
         pw1=2.0_dp*(lpd-m1);elw=lpd-pw
         mx=maxval(log_lik(:,j));w=exp(mx-log_lik(:,j));wmean=sum(w)/real(s,dp);w=w/wmean;w=min(w,sqrt(real(s,dp)))
         weighted=sum(exp(log_lik(:,j)-mx)*w);wsum=sum(w);eloo=mx+log(weighted/wsum);ploo=lpd-eloo
         res%pointwise(j,:)=[-2.0_dp*elw,-2.0_dp*(lpd-pw1),lpd,pw,pw1,elw,ploo,eloo]
      end do
      res%total=sum(res%pointwise,dim=1)
      do i=1,8
         m2=sum(res%pointwise(:,i))/real(n,dp);var=sum((res%pointwise(:,i)-m2)**2)/real(max(1,n-1),dp)
         res%se(i)=sqrt(real(n,dp)*var)
      end do
      res%waic=res%total(1);res%waic2=res%total(2);res%elpd_waic=res%total(6);res%p_waic=res%total(4);res%p_waic1=res%total(5)
      res%p_loo=res%total(7);res%elpd_loo=res%total(8);res%waic_ci=[res%waic-1.96_dp*res%se(1),res%waic+1.96_dp*res%se(1)]
   end function waic
end module mcmcpack_utils
