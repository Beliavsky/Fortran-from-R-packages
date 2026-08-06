! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran computational translation of tvgarch 2.4.3.
module tvgarch_tests
   use garchx_kinds, only : dp
   use garchx_linalg, only : solve_linear
   use tvgarch_model, only : tvgarch_spec, tvgarch_fit, make_tvgarch_spec, fit_tvgarch
   implicit none
   private

   type, public :: tvgarch_test_result
      type(tvgarch_fit) :: garch11
      real(dp) :: nonrobust(4,2) = 0.0_dp
      real(dp) :: robust(4,2) = 0.0_dp
      integer :: selected_order = 0
      real(dp), allocatable :: xtv(:)
      integer :: status = 1
   end type tvgarch_test_result
   public :: tvgarch_test
contains
   recursive pure real(dp) function log_gamma_lanczos(x) result(value)
      real(dp), intent(in) :: x
      real(dp), parameter :: c(9)=[0.99999999999980993_dp,676.5203681218851_dp,&
        -1259.1392167224028_dp,771.32342877765313_dp,-176.61502916214059_dp,&
        12.507343278686905_dp,-0.13857109526572012_dp,9.9843695780195716e-6_dp,&
        1.5056327351493116e-7_dp]
      real(dp)::z,t,s
      integer::i
      if(x<0.5_dp)then
         value=log(acos(-1.0_dp))-log(sin(acos(-1.0_dp)*x))-log_gamma_lanczos(1.0_dp-x)
      else
         z=x-1.0_dp;s=c(1)
         do i=2,size(c);s=s+c(i)/(z+real(i-1,dp));end do
         t=z+7.5_dp
         value=0.5_dp*log(2.0_dp*acos(-1.0_dp))+(z+0.5_dp)*log(t)-t+log(s)
      end if
   end function log_gamma_lanczos

   real(dp) function gamma_q(a,x) result(q)
      real(dp),intent(in)::a,x
      integer,parameter::itmax=200
      real(dp),parameter::eps=3.0e-14_dp,fpmin=1.0e-300_dp
      integer::n
      real(dp)::sumv,del,ap,b,c,d,h,an,p
      if(x<=0.0_dp)then;q=1.0_dp;return;end if
      if(x<a+1.0_dp)then
         ap=a;sumv=1.0_dp/a;del=sumv
         do n=1,itmax
            ap=ap+1.0_dp;del=del*x/ap;sumv=sumv+del
            if(abs(del)<abs(sumv)*eps)exit
         end do
         p=sumv*exp(-x+a*log(x)-log_gamma_lanczos(a));q=max(0.0_dp,1.0_dp-p)
      else
         b=x+1.0_dp-a;c=1.0_dp/fpmin;d=1.0_dp/b;h=d
         do n=1,itmax
            an=-real(n,dp)*(real(n,dp)-a);b=b+2.0_dp;d=an*d+b
            if(abs(d)<fpmin)d=fpmin;c=b+an/c;if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d;del=d*c;h=h*del;if(abs(del-1.0_dp)<eps)exit
         end do
         q=exp(-x+a*log(x)-log_gamma_lanczos(a))*h
      end if
   end function gamma_q

   real(dp) function chisq_upper(stat,df) result(p)
      real(dp),intent(in)::stat
      integer,intent(in)::df
      if(stat<=0.0_dp)then;p=1.0_dp;else;p=gamma_q(0.5_dp*real(df,dp),0.5_dp*stat);end if
   end function chisq_upper

   subroutine least_squares_residual(y,x,residual,status)
      real(dp),intent(in)::y(:,:),x(:,:)
      real(dp),allocatable,intent(out)::residual(:,:)
      integer,intent(out)::status
      real(dp),allocatable::xtx(:,:),xty(:,:),coef(:),rhs(:)
      integer::j,st
      if(size(y,1)/=size(x,1))then;status=1;allocate(residual(0,0));return;end if
      allocate(xtx(size(x,2),size(x,2)),xty(size(x,2),size(y,2)),residual(size(y,1),size(y,2)))
      xtx=matmul(transpose(x),x);xty=matmul(transpose(x),y)
      do j=1,size(y,2)
         rhs=xty(:,j);call solve_linear(xtx,rhs,coef,st)
         if(st/=0)then;status=2;residual=0.0_dp;return;end if
         residual(:,j)=y(:,j)-matmul(x,coef)
      end do
      status=0
   end subroutine least_squares_residual

   real(dp) function regression_rss(y,x,status) result(rss)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(out)::status
      real(dp),allocatable::ym(:,:),res(:,:)
      allocate(ym(size(y),1));ym(:,1)=y
      call least_squares_residual(ym,x,res,status)
      if(status/=0)then;rss=huge(1.0_dp);else;rss=sum(res(:,1)**2);end if
   end function regression_rss

   subroutine robust_lm(z,yaux,xbase,lm,pval,status)
      real(dp),intent(in)::z(:),yaux(:,:),xbase(:,:)
      real(dp),intent(out)::lm,pval
      integer,intent(out)::status
      real(dp),allocatable::res(:,:),weighted(:,:),ones(:),onesm(:,:),finalres(:,:)
      call least_squares_residual(yaux,xbase,res,status);if(status/=0)return
      allocate(weighted(size(res,1),size(res,2)),ones(size(z)),onesm(size(z),1))
      weighted=spread(z,2,size(res,2))*res;ones=1.0_dp;onesm(:,1)=ones
      call least_squares_residual(onesm,weighted,finalres,status);if(status/=0)return
      lm=real(size(z),dp)-sum(finalres(:,1)**2)
      pval=chisq_upper(lm,size(yaux,2))
   end subroutine robust_lm

   subroutine tvgarch_test(y,result,xtv,alpha,turbo)
      real(dp),intent(in)::y(:)
      type(tvgarch_test_result),intent(out)::result
      real(dp),intent(in),optional::xtv(:),alpha
      logical,intent(in),optional::turbo
      type(tvgarch_spec)::spec
      integer::n,i,st
      real(dp)::level,rss0,rss11,rss12,rss13,lm,p
      real(dp),allocatable::xv(:),y2(:),z(:),v(:,:),dhdw(:,:),base(:,:),x1(:,:),x2(:,:),x3(:,:)
      real(dp),allocatable::one(:),ya(:,:),xb(:,:)
      n=size(y);level=0.05_dp;if(present(alpha))level=alpha
      call make_tvgarch_spec(spec,order_h=[1,1,0],status=st)
      call fit_tvgarch(y,spec,result%garch11,turbo=turbo)
      if(.not.allocated(result%garch11%sigma2))then;result%status=1;return;end if
      allocate(xv(n));if(present(xtv))then
         if(size(xtv)/=n)then;result%status=2;return;end if;xv=xtv
      else
         do i=1,n;xv(i)=real(i,dp)/real(n,dp);end do
      end if
      allocate(result%xtv(n));result%xtv=xv
      allocate(y2(n),z(n),v(n,3),dhdw(n,3),base(n,3),one(n))
      y2=y*y;z=y2/result%garch11%sigma2-1.0_dp;one=1.0_dp
      v(1,:)=[1.0_dp,sum(y2)/real(n,dp),sum(y2)/real(n,dp)]
      do i=2,n;v(i,:)=[1.0_dp,y2(i-1),result%garch11%h(i-1)];end do
      dhdw(1,:)=v(1,:)
      do i=2,n;dhdw(i,:)=v(i,:)+result%garch11%hfit%par(3)*dhdw(i-1,:);end do
      base=spread(1.0_dp/result%garch11%h,2,3)*dhdw
      allocate(x1(n,5),x2(n,6),x3(n,7))
      x1(:,1:3)=base;x1(:,4)=one;x1(:,5)=xv
      x2(:,1:5)=x1;x2(:,6)=xv*xv
      x3(:,1:6)=x2;x3(:,7)=xv*xv*xv
      rss0=sum(z*z);rss11=regression_rss(z,x1,st);rss12=regression_rss(z,x2,st);rss13=regression_rss(z,x3,st)
      result%nonrobust(1,:)=[real(n,dp)*(rss0-rss13)/rss0,0.0_dp]
      result%nonrobust(2,:)=[real(n,dp)*(rss12-rss13)/rss12,0.0_dp]
      result%nonrobust(3,:)=[real(n,dp)*(rss11-rss12)/rss11,0.0_dp]
      result%nonrobust(4,:)=[real(n,dp)*(rss0-rss11)/rss0,0.0_dp]
      result%nonrobust(:,2)=[chisq_upper(result%nonrobust(1,1),3), &
         chisq_upper(result%nonrobust(2,1),1),chisq_upper(result%nonrobust(3,1),1), &
         chisq_upper(result%nonrobust(4,1),1)]
      allocate(ya(n,4));ya(:,1)=one;ya(:,2)=xv;ya(:,3)=xv*xv;ya(:,4)=xv*xv*xv
      call robust_lm(z,ya,base,lm,p,st);result%robust(1,:)=[lm,p]
      allocate(xb(n,6));xb(:,1:3)=base;xb(:,4)=one;xb(:,5)=xv;xb(:,6)=xv*xv
      ya=reshape([xv*xv*xv],[n,1]);call robust_lm(z,ya,xb,lm,p,st);result%robust(2,:)=[lm,p]
      deallocate(ya,xb);allocate(ya(n,1),xb(n,5));ya(:,1)=xv*xv
      xb(:,1:3)=base;xb(:,4)=one;xb(:,5)=xv
      call robust_lm(z,ya,xb,lm,p,st);result%robust(3,:)=[lm,p]
      deallocate(ya,xb);allocate(ya(n,2));ya(:,1)=one;ya(:,2)=xv
      call robust_lm(z,ya,base,lm,p,st);result%robust(4,:)=[lm,p]
      result%selected_order=0
      if(result%robust(1,2)<level)then
         if(result%robust(3,2)<result%robust(4,2) .and. result%robust(3,2)<result%robust(2,2) .and. &
            result%robust(3,2)<level)result%selected_order=2
         if(result%robust(4,2)<result%robust(2,2) .and. result%robust(4,2)<level)result%selected_order=1
         if(result%robust(2,2)<result%robust(4,2) .and. result%robust(2,2)<level)result%selected_order=3
      end if
      result%status=0
   end subroutine tvgarch_test
end module tvgarch_tests
