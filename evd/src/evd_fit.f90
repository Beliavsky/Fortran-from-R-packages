! SPDX-License-Identifier: GPL-3.0-only
module evd_fit
   use r_compat, only : dp, optim_result_t, optim_bfgs
   use evd_univariate, only : dgev, dgpd, dgumbelx
   use evd_bivariate, only : dbvlog, dbvalog, dbvhr, dbvneglog, dbvaneglog, &
      dbvbilog, dbvnegbilog, dbvct, dbvamix
   implicit none
   private
   public :: evd_fit_t, fit_gev, fit_gpd, fit_pp, fit_gumbelx, fit_bvevd

   type :: evd_fit_t
      real(dp), allocatable :: estimate(:)
      real(dp) :: nll = huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      integer :: convergence = 1
      integer :: counts(2) = 0
      real(dp), allocatable :: hessian(:,:)
      character(len=:), allocatable :: message
   end type evd_fit_t
contains

function fit_gev(x, start, maxit) result(fit)
   real(dp), intent(in) :: x(:)
   real(dp), intent(in), optional :: start(3)
   integer, intent(in), optional :: maxit
   type(evd_fit_t) :: fit
   type(optim_result_t) :: opt
   real(dp) :: p0(3), sc, mu, theta(3)
   integer :: mit
   sc=sqrt(max(6.0_dp*variance(x),epsilon(1.0_dp)))/acos(-1.0_dp)
   mu=sum(x)/real(size(x),dp)-0.58_dp*sc
   p0=[mu,sc,0.0_dp]
   if(present(start)) p0=start
   theta=[p0(1),log(max(p0(2),sqrt(tiny(1.0_dp)))),p0(3)]
   mit=300
   if(present(maxit)) mit=maxit
   opt=optim_bfgs(obj,theta,maxit=mit,reltol=1.0e-8_dp,ndeps=1.0e-5_dp,hessian=.true.)
   allocate(fit%estimate(3))
   fit%estimate=[opt%par(1),exp(opt%par(2)),opt%par(3)]
   call fill_fit(fit,opt)
contains
   pure function obj(p) result(v)
      real(dp),intent(in)::p(:)
      real(dp)::v,scale
      integer::i
      scale=exp(p(2))
      v=0.0_dp
      do i=1,size(x)
         v=v-dgev(x(i),p(1),scale,p(3),log_=.true.)
         if(.not.(v<1.0e100_dp)) then
         v=1.0e100_dp
         return
         end if
      end do
   end function obj
end function fit_gev

function fit_gpd(x, threshold, start, maxit) result(fit)
   real(dp), intent(in) :: x(:), threshold
   real(dp), intent(in), optional :: start(2)
   integer, intent(in), optional :: maxit
   type(evd_fit_t) :: fit
   type(optim_result_t) :: opt
   real(dp), allocatable :: exc(:)
   real(dp) :: p0(2), theta(2)
   integer :: i,n,mit
   n=count(x>threshold)
   allocate(exc(n))
   n=0
   do i=1,size(x)
   if(x(i)>threshold) then
   n=n+1
   exc(n)=x(i)
   end if
   end do
   if(n==0) then
      allocate(fit%estimate(2))
      fit%estimate=0.0_dp
      fit%message='no data above threshold'
      return
   end if
   p0=[max(sum(exc-threshold)/real(n,dp),sqrt(tiny(1.0_dp))),0.0_dp]
   if(present(start)) p0=start
   theta=[log(max(p0(1),sqrt(tiny(1.0_dp)))),p0(2)]
   mit=300
   if(present(maxit)) mit=maxit
   opt=optim_bfgs(obj,theta,maxit=mit,reltol=1.0e-8_dp,ndeps=1.0e-5_dp,hessian=.true.)
   allocate(fit%estimate(2))
   fit%estimate=[exp(opt%par(1)),opt%par(2)]
   call fill_fit(fit,opt)
contains
   pure function obj(p) result(v)
      real(dp),intent(in)::p(:)
      real(dp)::v,scale
      integer::j
      scale=exp(p(1))
      v=0.0_dp
      do j=1,size(exc)
         v=v-dgpd(exc(j),threshold,scale,p(2),log_=.true.)
         if(.not.(v<1.0e100_dp)) then
         v=1.0e100_dp
         return
         end if
      end do
   end function obj
end function fit_gpd

function fit_pp(x, threshold, npp, start, maxit) result(fit)
   real(dp), intent(in) :: x(:), threshold
   real(dp), intent(in), optional :: npp
   real(dp), intent(in), optional :: start(3)
   integer, intent(in), optional :: maxit
   type(evd_fit_t) :: fit
   type(optim_result_t) :: opt
   real(dp), allocatable :: exc(:)
   real(dp) :: p0(3), theta(3), nop, sc, mu
   integer :: i,n,mit
   n=count(x>threshold)
   allocate(exc(n))
   n=0
   do i=1,size(x)
   if(x(i)>threshold) then
   n=n+1
   exc(n)=x(i)
   end if
   end do
   if(n==0) then
      allocate(fit%estimate(3))
      fit%estimate=0.0_dp
      fit%message='no data above threshold'
      return
   end if
   nop=real(size(x),dp)
   if(present(npp)) nop=real(size(x),dp)/npp
   sc=sqrt(max(6.0_dp*variance(x),epsilon(1.0_dp)))/acos(-1.0_dp)
   mu=sum(x)/real(size(x),dp)+(log(nop)-0.58_dp)*sc
   p0=[mu,sc,0.0_dp]
   if(present(start)) p0=start
   theta=[p0(1),log(max(p0(2),sqrt(tiny(1.0_dp)))),p0(3)]
   mit=300
   if(present(maxit)) mit=maxit
   opt=optim_bfgs(obj,theta,maxit=mit,reltol=1.0e-8_dp,ndeps=1.0e-5_dp,hessian=.true.)
   allocate(fit%estimate(3))
   fit%estimate=[opt%par(1),exp(opt%par(2)),opt%par(3)]
   call fill_fit(fit,opt)
contains
   pure function obj(p) result(v)
      real(dp),intent(in)::p(:)
      real(dp)::v,scale,z,u,sh
      integer::j
      scale=exp(p(2))
      sh=p(3)
      v=0.0_dp
      do j=1,size(exc)
         z=(exc(j)-p(1))/scale
         if(abs(sh)<=epsilon(1.0_dp)**0.3_dp) then
            v=v+log(scale)+z
         else
            u=1.0_dp+sh*z
            if(u<=0.0_dp) then
            v=1.0e100_dp
            return
            end if
            v=v+log(scale)+(1.0_dp/sh+1.0_dp)*log(u)
         end if
      end do
      z=(threshold-p(1))/scale
      if(abs(sh)<=epsilon(1.0_dp)**0.3_dp) then
         v=v+nop*exp(-z)
      else
         u=1.0_dp+sh*z
         if(u<=0.0_dp) then
            if(sh>0.0_dp) v=1.0e100_dp
         else
            v=v+nop*u**(-1.0_dp/sh)
         end if
      end if
   end function obj
end function fit_pp

function fit_gumbelx(x, start, maxit) result(fit)
   real(dp), intent(in) :: x(:)
   real(dp), intent(in), optional :: start(4)
   integer, intent(in), optional :: maxit
   type(evd_fit_t) :: fit
   type(optim_result_t) :: opt
   real(dp) :: p0(4), theta(4), sc, mu
   integer :: mit
   sc=sqrt(max(6.0_dp*variance(x),epsilon(1.0_dp)))/acos(-1.0_dp)
   mu=sum(x)/real(size(x),dp)-0.58_dp*sc
   p0=[mu-sc/2.0_dp,sc,mu+sc/2.0_dp,sc]
   if(present(start)) p0=start
   theta=[p0(1),log(max(p0(2),sqrt(tiny(1.0_dp)))),log(max(p0(3)-p0(1),sqrt(tiny(1.0_dp)))), &
          log(max(p0(4),sqrt(tiny(1.0_dp))))]
   mit=400
   if(present(maxit)) mit=maxit
   opt=optim_bfgs(obj,theta,maxit=mit,reltol=1.0e-8_dp,ndeps=1.0e-5_dp,hessian=.true.)
   allocate(fit%estimate(4))
   fit%estimate=[opt%par(1),exp(opt%par(2)),opt%par(1)+exp(opt%par(3)),exp(opt%par(4))]
   call fill_fit(fit,opt)
contains
   pure function obj(p) result(v)
      real(dp),intent(in)::p(:)
      real(dp)::v,l1,s1,l2,s2
      integer::j
      l1=p(1)
      s1=exp(p(2))
      l2=l1+exp(p(3))
      s2=exp(p(4))
      v=0.0_dp
      do j=1,size(x)
         v=v-dgumbelx(x(j),l1,s1,l2,s2,log_=.true.)
         if(.not.(v<1.0e100_dp)) then
         v=1.0e100_dp
         return
         end if
      end do
   end function obj
end function fit_gumbelx

function fit_bvevd(data, model, start, maxit) result(fit)
   real(dp), intent(in) :: data(:,:)
   character(len=*), intent(in) :: model
   real(dp), intent(in) :: start(:)
   integer, intent(in), optional :: maxit
   type(evd_fit_t) :: fit
   type(optim_result_t) :: opt
   real(dp), allocatable :: theta(:)
   integer :: mit, m
   m=size(start)
   allocate(theta(m))
   theta=start
   call natural_to_theta(theta,model)
   mit=500
   if(present(maxit)) mit=maxit
   opt=optim_bfgs(obj,theta,maxit=mit,reltol=1.0e-8_dp,ndeps=2.0e-5_dp,hessian=.true.)
   allocate(fit%estimate(m))
   fit%estimate=opt%par
   call theta_to_natural(fit%estimate,model)
   call fill_fit(fit,opt)
contains
   pure function obj(p) result(v)
      real(dp),intent(in)::p(:)
      real(dp)::v, par(size(p)), ld
      integer::j
      par=p
      call theta_to_natural(par,model)
      v=0.0_dp
      do j=1,size(data,1)
         ld=biv_logdens(data(j,1),data(j,2),par,model)
         if(.not.(ld>-1.0e99_dp .and. ld<1.0e99_dp)) then
         v=1.0e100_dp
         return
         end if
         v=v-ld
      end do
   end function obj
end function fit_bvevd

pure function biv_logdens(x,y,p,model) result(v)
   real(dp),intent(in)::x,y,p(:)
   character(len=*),intent(in)::model
   real(dp)::v
   real(dp) :: mar1(3), mar2(3), asy(2)
   select case(trim(model))
   case('log')
      if(size(p)/=7) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(2:4)
      mar2=p(5:7)
      v=dbvlog(x,y,p(1),mar1,mar2,log_=.true.)
   case('hr')
      if(size(p)/=7) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(2:4)
      mar2=p(5:7)
      v=dbvhr(x,y,p(1),mar1,mar2,log_=.true.)
   case('neglog')
      if(size(p)/=7) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(2:4)
      mar2=p(5:7)
      v=dbvneglog(x,y,p(1),mar1,mar2,log_=.true.)
   case('alog')
      if(size(p)/=9) then
      v=-huge(1.0_dp)
      return
      end if
      asy=p(2:3)
      mar1=p(4:6)
      mar2=p(7:9)
      v=dbvalog(x,y,p(1),asy,mar1,mar2,log_=.true.)
   case('aneglog')
      if(size(p)/=9) then
      v=-huge(1.0_dp)
      return
      end if
      asy=p(2:3)
      mar1=p(4:6)
      mar2=p(7:9)
      v=dbvaneglog(x,y,p(1),asy,mar1,mar2,log_=.true.)
   case('bilog')
      if(size(p)/=8) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(3:5)
      mar2=p(6:8)
      v=dbvbilog(x,y,p(1),p(2),mar1,mar2,log_=.true.)
   case('negbilog')
      if(size(p)/=8) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(3:5)
      mar2=p(6:8)
      v=dbvnegbilog(x,y,p(1),p(2),mar1,mar2,log_=.true.)
   case('ct')
      if(size(p)/=8) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(3:5)
      mar2=p(6:8)
      v=dbvct(x,y,p(1),p(2),mar1,mar2,log_=.true.)
   case('amix')
      if(size(p)/=8) then
      v=-huge(1.0_dp)
      return
      end if
      mar1=p(3:5)
      mar2=p(6:8)
      v=dbvamix(x,y,p(1),p(2),mar1,mar2,log_=.true.)
   case default
      v=-huge(1.0_dp)
   end select
end function biv_logdens

pure subroutine natural_to_theta(p,model)
   real(dp),intent(inout)::p(:)
   character(len=*),intent(in)::model
   select case(trim(model))
   case('log')
   p(1)=logit(p(1))
   call margins_to_theta(p,2)
   case('hr','neglog')
   p(1)=log(max(p(1),tiny(1.0_dp)))
   call margins_to_theta(p,2)
   case('alog')
   p(1)=logit(p(1))
   p(2)=logit(p(2))
   p(3)=logit(p(3))
   call margins_to_theta(p,4)
   case('aneglog')
   p(1)=log(max(p(1),tiny(1.0_dp)))
   p(2)=logit(p(2))
   p(3)=logit(p(3))
   call margins_to_theta(p,4)
   case('bilog')
   p(1)=logit(p(1))
   p(2)=logit(p(2))
   call margins_to_theta(p,3)
   case('negbilog','ct')
   p(1)=log(max(p(1),tiny(1.0_dp)))
   p(2)=log(max(p(2),tiny(1.0_dp)))
   call margins_to_theta(p,3)
   case('amix'); call margins_to_theta(p,3)
   end select
end subroutine natural_to_theta

pure subroutine theta_to_natural(p,model)
   real(dp),intent(inout)::p(:)
   character(len=*),intent(in)::model
   select case(trim(model))
   case('log')
   p(1)=logistic(p(1))
   call margins_from_theta(p,2)
   case('hr','neglog')
   p(1)=exp(p(1))
   call margins_from_theta(p,2)
   case('alog')
   p(1)=logistic(p(1))
   p(2)=logistic(p(2))
   p(3)=logistic(p(3))
   call margins_from_theta(p,4)
   case('aneglog')
   p(1)=exp(p(1))
   p(2)=logistic(p(2))
   p(3)=logistic(p(3))
   call margins_from_theta(p,4)
   case('bilog')
   p(1)=logistic(p(1))
   p(2)=logistic(p(2))
   call margins_from_theta(p,3)
   case('negbilog','ct')
   p(1)=exp(p(1))
   p(2)=exp(p(2))
   call margins_from_theta(p,3)
   case('amix'); call margins_from_theta(p,3)
   end select
end subroutine theta_to_natural

pure subroutine margins_to_theta(p,k)
   real(dp),intent(inout)::p(:)
   integer,intent(in)::k
   if(size(p)>=k+5) then
   p(k+1)=log(max(p(k+1),tiny(1.0_dp)))
   p(k+4)=log(max(p(k+4),tiny(1.0_dp)))
   end if
end subroutine margins_to_theta
pure subroutine margins_from_theta(p,k)
   real(dp),intent(inout)::p(:)
   integer,intent(in)::k
   if(size(p)>=k+5) then
   p(k+1)=exp(p(k+1))
   p(k+4)=exp(p(k+4))
   end if
end subroutine margins_from_theta

pure elemental function logistic(x) result(v)
   real(dp),intent(in)::x
   real(dp)::v
   if(x>=0.0_dp) then
   v=1.0_dp/(1.0_dp+exp(-x))
   else
   v=exp(x)/(1.0_dp+exp(x))
   end if
end function logistic
pure elemental function logit(x) result(v)
   real(dp),intent(in)::x
   real(dp)::v
   v=log(max(x,tiny(1.0_dp))/max(1.0_dp-x,tiny(1.0_dp)))
end function logit

pure function variance(x) result(v)
   real(dp),intent(in)::x(:)
   real(dp)::v,m
   if(size(x)<2) then
   v=0.0_dp
   return
   end if
   m=sum(x)/real(size(x),dp)
   v=sum((x-m)**2)/real(size(x)-1,dp)
end function variance

subroutine fill_fit(fit,opt)
   type(evd_fit_t),intent(inout)::fit
   type(optim_result_t),intent(in)::opt
   fit%nll=opt%value
   fit%deviance=2.0_dp*opt%value
   fit%convergence=opt%convergence
   fit%counts=opt%counts
   if(allocated(opt%hessian)) then
   allocate(fit%hessian(size(opt%hessian,1),size(opt%hessian,2)))
   fit%hessian=opt%hessian
   end if
   if(allocated(opt%message)) fit%message=opt%message
end subroutine fill_fit
end module evd_fit
