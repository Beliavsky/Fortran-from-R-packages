! SPDX-License-Identifier: MIT
module gkwdist_startvalues
   use, intrinsic :: iso_fortran_env, only : int64
   use gkwdist_kinds, only : dp
   use gkwdist_math, only : finite_dp
   use gkwdist_core, only : fam_gkw,fam_bkw,fam_kkw,fam_ekw,fam_mc,fam_kw,fam_beta, &
      family_from_name,family_npar,family_full_parameters,dgkw_scalar
   implicit none
   private
   public :: gkwgetstartvalues, theoretical_moment
contains

   pure function family_pdf(x,theta,fam) result(f)
      real(dp),intent(in)::x,theta(:)
      integer,intent(in)::fam
      real(dp)::f,a,b,g,d,l
      logical::ok
      call family_full_parameters(fam,theta,a,b,g,d,l,ok)
      if(ok) then
         f=dgkw_scalar(x,a,b,g,d,l)
      else
         f=0.0_dp
      end if
   end function family_pdf

   pure function theoretical_moment(r,theta,fam) result(mom)
      integer,intent(in)::r,fam
      real(dp),intent(in)::theta(:)
      real(dp)::mom,h,x,fx,sumv
      integer,parameter::np=51
      integer::i,w
      h=1.0_dp/real(np-1,dp); sumv=0.0_dp
      do i=0,np-1
         x=real(i,dp)*h
         if(i==0 .or. i==np-1) then
            w=1
         else if(mod(i,2)==1) then
            w=4
         else
            w=2
         end if
         fx=x**r*family_pdf(x,theta,fam)
         sumv=sumv+real(w,dp)*fx
      end do
      mom=h*sumv/3.0_dp
      if(.not.finite_dp(mom) .or. abs(mom)<1.0e-14_dp) mom=0.5_dp
   end function theoretical_moment

   pure function objective(theta,sample_mom,fam) result(obj)
      real(dp),intent(in)::theta(:),sample_mom(5)
      integer,intent(in)::fam
      real(dp)::obj,theor,rel
      real(dp),parameter::weights(5)=[1.0_dp,0.8_dp,0.6_dp,0.4_dp,0.2_dp]
      integer::r
      if(any(theta<=0.0_dp) .or. any(.not.[(finite_dp(theta(r)),r=1,size(theta))])) then
         obj=huge(1.0_dp); return
      end if
      obj=0.0_dp
      do r=1,5
         theor=theoretical_moment(r,theta,fam)
         if(abs(sample_mom(r))<1.0e-10_dp) then
            obj=obj+weights(r)*theor*theor
         else
            rel=(theor-sample_mom(r))/sample_mom(r)
            obj=obj+weights(r)*rel*rel
         end if
      end do
      if(.not.finite_dp(obj)) obj=huge(1.0_dp)
   end function objective

   subroutine sort_simplex(simplex,vals)
      real(dp),intent(inout)::simplex(:,:),vals(:)
      integer::i,j,k,n
      real(dp)::tv
      real(dp),allocatable::col(:)
      n=size(vals); allocate(col(size(simplex,1)))
      do i=1,n-1
         k=i
         do j=i+1,n
            if(vals(j)<vals(k)) k=j
         end do
         if(k/=i) then
            tv=vals(i); vals(i)=vals(k); vals(k)=tv
            col=simplex(:,i); simplex(:,i)=simplex(:,k); simplex(:,k)=col
         end if
      end do
   end subroutine sort_simplex

   function nelder_mead(initial,sample_mom,fam,max_iter,tol) result(best)
      real(dp),intent(in)::initial(:),sample_mom(5)
      integer,intent(in)::fam
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol
      real(dp),allocatable::best(:)
      real(dp),allocatable::simplex(:,:),vals(:),centroid(:),reflected(:),expanded(:),contracted(:)
      real(dp)::fr,fe,fc,diameter,step,tolerance
      integer::n,i,it,niter
      n=size(initial); niter=1000; tolerance=1.0e-6_dp
      if(present(max_iter)) niter=max_iter
      if(present(tol)) tolerance=tol
      allocate(simplex(n,n+1),vals(n+1),centroid(n),reflected(n),expanded(n),contracted(n),best(n))
      simplex(:,1)=initial
      do i=2,n+1
         simplex(:,i)=initial
         step=0.05_dp*max(abs(initial(i-1)),0.1_dp)
         simplex(i-1,i)=simplex(i-1,i)+step
         where(simplex(:,i)<=0.0_dp) simplex(:,i)=0.01_dp
      end do
      do i=1,n+1; vals(i)=objective(simplex(:,i),sample_mom,fam); end do
      do it=1,niter
         call sort_simplex(simplex,vals)
         diameter=0.0_dp
         do i=2,n+1
            diameter=max(diameter,sqrt(sum((simplex(:,i)-simplex(:,1))**2)))
         end do
         if(diameter<tolerance) exit
         centroid=0.0_dp
         do i=1,n; centroid=centroid+simplex(:,i); end do
         centroid=centroid/real(n,dp)
         reflected=centroid+(centroid-simplex(:,n+1)); where(reflected<=0.0_dp) reflected=0.01_dp
         fr=objective(reflected,sample_mom,fam)
         if(fr<vals(1)) then
            expanded=centroid+2.0_dp*(reflected-centroid); where(expanded<=0.0_dp) expanded=0.01_dp
            fe=objective(expanded,sample_mom,fam)
            if(fe<fr) then; simplex(:,n+1)=expanded; vals(n+1)=fe
            else; simplex(:,n+1)=reflected; vals(n+1)=fr; end if
         else if(fr<vals(n)) then
            simplex(:,n+1)=reflected; vals(n+1)=fr
         else
            if(fr<vals(n+1)) then
               contracted=centroid+0.5_dp*(reflected-centroid)
               where(contracted<=0.0_dp) contracted=0.01_dp
               fc=objective(contracted,sample_mom,fam)
               if(fc<=fr) then
                  simplex(:,n+1)=contracted; vals(n+1)=fc
                  cycle
               end if
            else
               contracted=centroid-0.5_dp*(centroid-simplex(:,n+1))
               where(contracted<=0.0_dp) contracted=0.01_dp
               fc=objective(contracted,sample_mom,fam)
               if(fc<vals(n+1)) then
                  simplex(:,n+1)=contracted; vals(n+1)=fc
                  cycle
               end if
            end if
            do i=2,n+1
               simplex(:,i)=simplex(:,1)+0.5_dp*(simplex(:,i)-simplex(:,1))
               where(simplex(:,i)<=0.0_dp) simplex(:,i)=0.01_dp
               vals(i)=objective(simplex(:,i),sample_mom,fam)
            end do
         end if
      end do
      call sort_simplex(simplex,vals); best=simplex(:,1)
   end function nelder_mead

   function lcg_uniform(state) result(u)
      integer(int64),intent(inout)::state
      real(dp)::u
      integer(int64)::hi,lo,test
      hi=state/127773_int64; lo=mod(state,127773_int64)
      test=16807_int64*lo-2836_int64*hi
      if(test<=0_int64) test=test+2147483647_int64
      state=test; u=real(state,dp)/2147483647.0_dp
   end function lcg_uniform

   subroutine make_initial(start,fam,k,alpha0,beta0,m1,state)
      real(dp),intent(out)::start(:)
      integer,intent(in)::fam,k
      real(dp),intent(in)::alpha0,beta0,m1
      integer(int64),intent(inout)::state
      real(dp)::u
      select case(fam)
      case(fam_gkw)
         select case(k)
         case(1); start=[alpha0,beta0,1.0_dp,0.1_dp,1.0_dp]
         case(2); start=[2.0_dp,2.0_dp,1.0_dp,0.5_dp,1.0_dp]
         case(3); start=[1.0_dp,1.0_dp,1.0_dp,0.1_dp,1.0_dp]
         case(4); start=[4.0_dp,2.0_dp,0.8_dp,0.5_dp,1.0_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(3)=0.5_dp+1.5_dp*u
            u=lcg_uniform(state); start(4)=0.1_dp+0.9_dp*u
            u=lcg_uniform(state); start(5)=0.5_dp+1.5_dp*u
         end select
      case(fam_bkw)
         select case(k)
         case(1); start=[alpha0,beta0,1.0_dp,0.5_dp]
         case(2); start=[2.0_dp,2.0_dp,1.0_dp,0.5_dp]
         case(3); start=[1.0_dp,1.0_dp,0.8_dp,0.3_dp]
         case(4); start=[3.0_dp,2.0_dp,1.5_dp,0.5_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(3)=0.5_dp+1.5_dp*u
            u=lcg_uniform(state); start(4)=0.1_dp+0.9_dp*u
         end select
      case(fam_kkw)
         select case(k)
         case(1); start=[alpha0,beta0,0.5_dp,1.0_dp]
         case(2); start=[2.0_dp,2.0_dp,0.5_dp,1.0_dp]
         case(3); start=[1.0_dp,1.0_dp,0.3_dp,1.2_dp]
         case(4); start=[3.0_dp,2.0_dp,0.7_dp,1.5_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(3)=0.1_dp+0.9_dp*u
            u=lcg_uniform(state); start(4)=0.5_dp+1.5_dp*u
         end select
      case(fam_ekw)
         select case(k)
         case(1); start=[alpha0,beta0,1.0_dp]
         case(2); start=[2.0_dp,2.0_dp,1.0_dp]
         case(3); start=[1.0_dp,1.0_dp,1.2_dp]
         case(4); start=[3.0_dp,2.0_dp,1.5_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(3)=0.5_dp+1.5_dp*u
         end select
      case(fam_mc)
         select case(k)
         case(1); start=[merge(2.0_dp,1.0_dp,m1>0.5_dp),merge(2.0_dp,1.0_dp,m1<0.5_dp),1.0_dp]
         case(2); start=[1.0_dp,1.0_dp,1.0_dp]
         case(3); start=[2.0_dp,2.0_dp,1.2_dp]
         case(4); start=[1.5_dp,1.5_dp,1.5_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+4.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+4.5_dp*u
            u=lcg_uniform(state); start(3)=0.5_dp+1.5_dp*u
         end select
      case(fam_kw)
         select case(k)
         case(1); start=[alpha0,beta0]
         case(2); start=[2.0_dp,2.0_dp]
         case(3); start=[1.0_dp,1.0_dp]
         case(4); start=[3.0_dp,2.0_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
         end select
      case(fam_beta)
         select case(k)
         case(1); start=[alpha0,beta0]
         case(2); start=[2.0_dp,2.0_dp]
         case(3); start=[1.0_dp,1.0_dp]
         case(4); start=[3.0_dp,2.0_dp]
         case default
            u=lcg_uniform(state); start(1)=0.5_dp+9.5_dp*u
            u=lcg_uniform(state); start(2)=0.5_dp+9.5_dp*u
         end select
      end select
   end subroutine make_initial

   subroutine constrain(theta,fam)
      real(dp),intent(inout)::theta(:)
      integer,intent(in)::fam
      integer :: k
      where(.not.[(finite_dp(theta(k)),k=1,size(theta))] .or. theta<=0.0_dp) theta=1.0_dp
      select case(fam)
      case(fam_gkw)
         theta(1)=min(50.0_dp,max(0.1_dp,theta(1))); theta(2)=min(50.0_dp,max(0.1_dp,theta(2)))
         theta(3)=min(10.0_dp,max(0.1_dp,theta(3))); theta(4)=min(10.0_dp,max(0.01_dp,theta(4)))
         theta(5)=min(20.0_dp,max(0.1_dp,theta(5)))
      case(fam_bkw)
         theta(1)=min(50.0_dp,max(0.1_dp,theta(1))); theta(2)=min(50.0_dp,max(0.1_dp,theta(2)))
         theta(3)=min(10.0_dp,max(0.1_dp,theta(3))); theta(4)=min(10.0_dp,max(0.01_dp,theta(4)))
      case(fam_kkw)
         theta(1)=min(50.0_dp,max(0.1_dp,theta(1))); theta(2)=min(50.0_dp,max(0.1_dp,theta(2)))
         theta(3)=min(10.0_dp,max(0.01_dp,theta(3))); theta(4)=min(20.0_dp,max(0.1_dp,theta(4)))
      case(fam_ekw)
         theta(1)=min(50.0_dp,max(0.1_dp,theta(1))); theta(2)=min(50.0_dp,max(0.1_dp,theta(2)))
         theta(3)=min(20.0_dp,max(0.1_dp,theta(3)))
      case(fam_mc)
         theta=min(20.0_dp,max(0.1_dp,theta))
      case(fam_kw,fam_beta)
         theta=min(50.0_dp,max(0.1_dp,theta))
      end select
   end subroutine constrain

   function gkwgetstartvalues(x,family,n_starts,objective_value,success) result(par)
      real(dp),intent(in)::x(:)
      character(len=*),intent(in),optional::family
      integer,intent(in),optional::n_starts
      real(dp),intent(out),optional::objective_value
      logical,intent(out),optional::success
      real(dp),allocatable::par(:)
      real(dp),allocatable::clean(:),start(:),opt(:)
      real(dp)::sm(5),m1,m2,var,alpha0,beta0,obj,best_obj
      integer::fam,np,ns,k,nclean,i
      integer(int64)::state
      logical::found
      fam=fam_gkw; if(present(family)) fam=family_from_name(family)
      np=family_npar(fam)
      if(np<=0) then
         allocate(par(0)); if(present(success)) success=.false.; if(present(objective_value)) objective_value=huge(1.0_dp)
         return
      end if
      allocate(par(np)); par=1.0_dp
      nclean=count([(finite_dp(x(i)),i=1,size(x))])
      if(nclean==0) then
         if(present(success)) success=.false.; if(present(objective_value)) objective_value=huge(1.0_dp); return
      end if
      allocate(clean(nclean)); k=0
      do i=1,size(x)
         if(finite_dp(x(i))) then
            k=k+1; clean(k)=max(1.0e-10_dp,min(1.0_dp-1.0e-10_dp,x(i)))
         end if
      end do
      do k=1,5; sm(k)=sum(clean**k)/real(nclean,dp); end do
      m1=max(0.01_dp,min(0.99_dp,sm(1))); m2=sm(2); var=m2-m1*m1
      if(var<=1.0e-10_dp) var=0.01_dp
      alpha0=max(0.1_dp,m1*(1.0_dp-m1)/var-m1)
      beta0=max(0.1_dp,alpha0*m1/(1.0_dp-m1))
      alpha0=min(20.0_dp,max(0.1_dp,alpha0)); beta0=min(20.0_dp,max(0.1_dp,beta0))
      if(fam==fam_beta) then
         alpha0=m1*((m1*(1.0_dp-m1)/var)-1.0_dp)
         beta0=(1.0_dp-m1)*((m1*(1.0_dp-m1)/var)-1.0_dp)
         alpha0=min(50.0_dp,max(0.1_dp,alpha0)); beta0=min(50.0_dp,max(0.1_dp,beta0))
      end if
      ns=5; if(present(n_starts)) ns=max(1,n_starts); ns=max(4,ns)
      allocate(start(np)); state=42_int64; best_obj=huge(1.0_dp); found=.false.
      do k=1,ns
         call make_initial(start,fam,k,alpha0,beta0,m1,state)
         obj=objective(start,sm,fam)
         if(obj<best_obj) then
            opt=nelder_mead(start,sm,fam)
            obj=objective(opt,sm,fam)
            if(obj<best_obj .and. finite_dp(obj)) then
               par=opt; best_obj=obj; found=.true.
            end if
         end if
      end do
      if(.not.found) then
         select case(fam)
         case(fam_gkw); par=[1.0_dp,1.0_dp,1.0_dp,0.1_dp,1.0_dp]
         case(fam_bkw); par=[1.0_dp,1.0_dp,1.0_dp,0.5_dp]
         case(fam_kkw); par=[1.0_dp,1.0_dp,0.5_dp,1.0_dp]
         case(fam_ekw,fam_mc); par=[1.0_dp,1.0_dp,1.0_dp]
         case(fam_kw,fam_beta); par=[1.0_dp,1.0_dp]
         end select
         best_obj=objective(par,sm,fam)
      end if
      call constrain(par,fam)
      if(present(objective_value)) objective_value=best_obj
      if(present(success)) success=found
   end function gkwgetstartvalues

end module gkwdist_startvalues
