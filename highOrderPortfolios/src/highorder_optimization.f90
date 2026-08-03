! SPDX-License-Identifier: GPL-3.0-only
module highorder_optimization
   use fitheavytail_kinds, only: dp
   use highorder_types
   use highorder_linalg
   use highorder_moments, only: evaluate_sample_details, evaluate_skew_t_details
   implicit none
   private

   public :: design_mvsk_portfolio_via_sample_moments
   public :: design_mvsk_portfolio_via_skew_t
   public :: design_mvsktilting_portfolio_via_sample_moments

contains

   subroutine design_mvsk_portfolio_via_sample_moments(lmd,statistics,result, &
      w_init,leverage,method,tau_w,gamma,zeta,maxiter,ftol,wtol,stopval)
      real(dp), intent(in) :: lmd(4)
      type(sample_moments), intent(in) :: statistics
      type(portfolio_result), intent(out) :: result
      real(dp), intent(in), optional :: w_init(:),leverage,tau_w,gamma,zeta
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: ftol,wtol,stopval
      real(dp), allocatable :: w(:),wold(:),what(:),g(:,:),h3(:,:),h4(:,:)
      real(dp), allocatable :: q(:,:),qvec(:),hn(:,:),hpsd(:,:),hist(:)
      real(dp) :: m(4),mnew(4),obj,objnew,tauw,gam,zet,ft,wt,sv,lev,rho
      character(len=16) :: meth
      integer :: n,iter,niter,i
      logical :: wc,fc,sc

      call clear_portfolio_result(result)
      n=statistics%nassets
      if(n<1 .or. statistics%status/=hop_success) then
         call fail_result(result,hop_invalid_argument,'invalid sample moments')
         return
      end if
      lev=1.0_dp
      if(present(leverage)) lev=leverage
      if(abs(lev-1.0_dp)>1.0e-12_dp) then
         call fail_result(result,hop_invalid_argument,'only leverage=1 is supported')
         return
      end if
      meth='Q-MVSK'
      if(present(method)) meth=to_upper(trim(method))
      if(meth/='Q-MVSK' .and. meth/='MM' .and. meth/='DC') then
         call fail_result(result,hop_invalid_argument,'method must be Q-MVSK, MM, or DC')
         return
      end if
      tauw=0.0_dp
      if(present(tau_w)) tauw=max(0.0_dp,tau_w)
      gam=1.0_dp
      if(present(gamma)) gam=gamma
      zet=1.0e-8_dp
      if(present(zeta)) zet=zeta
      niter=100
      if(present(maxiter)) niter=maxiter
      ft=1.0e-5_dp
      if(present(ftol)) ft=ftol
      wt=1.0e-4_dp
      if(present(wtol)) wt=wtol
      sv=-huge(1.0_dp)
      if(present(stopval)) sv=stopval
      if(niter<1 .or. gam<=0.0_dp .or. gam>1.0_dp .or. zet<0.0_dp) then
         call fail_result(result,hop_invalid_argument,'invalid optimization controls')
         return
      end if

      allocate(w(n),wold(n),what(n),g(4,n),h3(n,n),h4(n,n))
      allocate(q(n,n),qvec(n),hn(n,n),hpsd(n,n),hist(niter+1))
      if(present(w_init)) then
         if(size(w_init)/=n) then
            call fail_result(result,hop_dimension_mismatch,'w_init has wrong length')
            return
         end if
         call project_simplex(w_init,w)
      else
         w=1.0_dp/real(n,dp)
      end if
      call evaluate_sample_details(w,statistics,m,g,h3,h4)
      obj=mvsk_objective(lmd,m)
      hist(1)=obj

      do iter=1,niter
         wold=w
         select case(trim(meth))
         case('Q-MVSK')
            hn=-lmd(3)*h3+lmd(4)*h4
            call psd_projection(hn,hpsd)
            q=2.0_dp*lmd(2)*statistics%covariance+hpsd
            do i=1,n
               q(i,i)=q(i,i)+tauw+1.0e-10_dp
            end do
            qvec=lmd(1)*g(1,:)+lmd(3)*g(3,:)-lmd(4)*g(4,:) + &
                 matmul(hpsd,w)+tauw*w
         case('MM')
            rho=spectral_radius(-lmd(3)*h3+lmd(4)*h4)+1.0e-8_dp
            q=2.0_dp*lmd(2)*statistics%covariance
            do i=1,n
               q(i,i)=q(i,i)+rho
            end do
            qvec=lmd(1)*g(1,:)+lmd(3)*g(3,:)-lmd(4)*g(4,:)+rho*w
         case('DC')
            rho=spectral_radius(2.0_dp*lmd(2)*statistics%covariance - &
                 lmd(3)*h3+lmd(4)*h4)+1.0e-8_dp
            q=0.0_dp
            do i=1,n
               q(i,i)=rho
            end do
            qvec=rho*w+lmd(1)*g(1,:)-lmd(2)*g(2,:) + &
                 lmd(3)*g(3,:)-lmd(4)*g(4,:)
         end select
         call solve_simplex_qp(q,qvec,w,what,max_iter=1000,tol=1.0e-11_dp)
         w=w+gam*(what-w)
         call project_simplex(w,what)
         w=what
         gam=gam*(1.0_dp-zet*gam)
         call evaluate_sample_details(w,statistics,mnew,g,h3,h4)
         objnew=mvsk_objective(lmd,mnew)
         hist(iter+1)=objnew
         wc=vector_norm2(w-wold)<=wt*max(1.0_dp,vector_norm2(wold))
         fc=abs(objnew-obj)<=0.5_dp*ft*(abs(objnew)+abs(obj)+1.0e-16_dp)
         sc=objnew<=sv
         obj=objnew
         m=mnew
         if(wc .or. fc .or. sc) exit
      end do
      call finish_result(result,w,m,hist,iter,obj,iter<niter)
   end subroutine design_mvsk_portfolio_via_sample_moments

   subroutine design_mvsk_portfolio_via_skew_t(lambda,parameters,result, &
      w_init,method,gamma,zeta,tau_w,beta,tau,initial_eta,maxiter,ftol,wtol,stopval)
      real(dp), intent(in) :: lambda(4)
      type(skew_t_parameters), intent(in) :: parameters
      type(portfolio_result), intent(out) :: result
      real(dp), intent(in), optional :: w_init(:),gamma,zeta,tau_w,beta,tau,initial_eta
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: ftol,wtol,stopval
      real(dp), allocatable :: w(:),wold(:),what(:),w1(:),w2(:),r(:),vvec(:)
      real(dp), allocatable :: g(:,:),g1(:,:),h2(:,:),h3(:,:),h4(:,:),q(:,:),qvec(:)
      real(dp), allocatable :: hn(:,:),hpsd(:,:),hist(:),grad(:)
      real(dp) :: lam(4),scale,m(4),mnew(4),obj,objnew,gam,zet,tauw,bet,tauv,eta
      real(dp) :: ft,wt,sv,rho,alpha,denom
      character(len=16) :: meth
      integer :: n,niter,iter,i
      logical :: wc,fc,sc

      call clear_portfolio_result(result)
      if((parameters%status/=hop_success .and. parameters%status/=hop_not_converged) .or. &
         .not.allocated(parameters%mu)) then
         call fail_result(result,hop_invalid_argument,'invalid skew-t parameters')
         return
      end if
      n=size(parameters%mu)
      scale=maxval(abs(lambda))
      if(scale<=0.0_dp) then
         call fail_result(result,hop_invalid_argument,'lambda must not be all zero')
         return
      end if
      lam=lambda/scale
      meth='L-MVSK'
      if(present(method)) meth=to_upper(trim(method))
      if(.not.valid_skew_method(meth)) then
         call fail_result(result,hop_invalid_argument,'unknown skew-t MVSK method')
         return
      end if
      gam=1.0_dp
      if(present(gamma)) gam=gamma
      zet=1.0e-8_dp
      if(present(zeta)) zet=zeta
      tauw=0.0_dp
      if(present(tau_w)) tauw=max(0.0_dp,tau_w)
      bet=0.5_dp
      if(present(beta)) bet=beta
      tauv=1.0e5_dp
      if(present(tau)) tauv=tau
      eta=5.0_dp
      if(present(initial_eta)) eta=initial_eta
      niter=1000
      if(present(maxiter)) niter=maxiter
      ft=1.0e-6_dp
      if(present(ftol)) ft=ftol
      wt=1.0e-6_dp
      if(present(wtol)) wt=wtol
      sv=-huge(1.0_dp)
      if(present(stopval)) sv=stopval
      if(niter<1 .or. bet<=0.0_dp .or. bet>=1.0_dp .or. tauv<=0.0_dp) then
         call fail_result(result,hop_invalid_argument,'invalid optimization controls')
         return
      end if

      allocate(w(n),wold(n),what(n),w1(n),w2(n),r(n),vvec(n),grad(n))
      allocate(g(4,n),g1(4,n),h2(n,n),h3(n,n),h4(n,n),q(n,n),qvec(n))
      allocate(hn(n,n),hpsd(n,n),hist(niter+1))
      if(present(w_init)) then
         if(size(w_init)/=n) then
            call fail_result(result,hop_dimension_mismatch,'w_init has wrong length')
            return
         end if
         call project_simplex(w_init,w)
      else
         w=1.0_dp/real(n,dp)
      end if
      call evaluate_skew_t_details(w,parameters,m,g,h2,h3,h4)
      obj=mvsk_objective(lam,m)
      hist(1)=obj

      do iter=1,niter
         wold=w
         select case(trim(meth))
         case('Q-MVSK')
            hn=-lam(3)*h3+lam(4)*h4
            call psd_projection(hn,hpsd)
            q=lam(2)*h2+hpsd
            do i=1,n
               q(i,i)=q(i,i)+tauw+1.0e-10_dp
            end do
            qvec=lam(1)*g(1,:)+lam(3)*g(3,:)-lam(4)*g(4,:) + &
                 matmul(hpsd,w)+tauw*w
            call solve_simplex_qp(q,qvec,w,what,max_iter=1000,tol=1.0e-11_dp)
            w=w+gam*(what-w)
            gam=gam*(1.0_dp-zet*gam)
         case('L-MVSK')
            rho=spectral_radius(-lam(3)*h3+lam(4)*h4)+1.0e-8_dp
            q=lam(2)*h2
            do i=1,n
               q(i,i)=q(i,i)+rho
            end do
            qvec=rho*w+lam(1)*g(1,:)+lam(3)*g(3,:)-lam(4)*g(4,:)
            call solve_simplex_qp(q,qvec,w,what,max_iter=1000,tol=1.0e-11_dp)
            w=what
         case('DC')
            rho=spectral_radius(lam(2)*h2-lam(3)*h3+lam(4)*h4)+1.0e-8_dp
            q=0.0_dp
            do i=1,n
               q(i,i)=rho
            end do
            qvec=rho*w+lam(1)*g(1,:)-lam(2)*g(2,:) + &
                 lam(3)*g(3,:)-lam(4)*g(4,:)
            call solve_simplex_qp(q,qvec,w,what,max_iter=1000,tol=1.0e-11_dp)
            w=what
         case('PGD')
            grad=objective_gradient(lam,g)
            call pgd_backtrack_skew(w,grad,obj,parameters,lam,eta,bet,what,objnew)
            w=what
         case('RFPA','SQUAREM')
            grad=objective_gradient(lam,g)
            call project_simplex(w-grad/tauv,w1)
            call evaluate_skew_t_details(w1,parameters,mnew,g1)
            call project_simplex(w1-objective_gradient(lam,g1)/tauv,w2)
            r=w1-w
            vvec=(w2-w1)-r
            denom=vector_norm2(vvec)
            if(denom<=1.0e-15_dp) then
               alpha=-1.0_dp
            else
               alpha=-vector_norm2(r)/denom
            end if
            if(trim(meth)=='SQUAREM') alpha=min(-1.0_dp,alpha)
            if(trim(meth)=='RFPA') then
               denom=dot_product(r,vvec)
               if(denom<0.0_dp) alpha=max(alpha,dot_product(r,r)/denom)
            end if
            what=w-2.0_dp*alpha*r+alpha*alpha*vvec
            call evaluate_skew_t_details(what,parameters,mnew,g1)
            call project_simplex(what-objective_gradient(lam,g1)/tauv,w1)
            what=w1
            call evaluate_skew_t_details(what,parameters,mnew)
            objnew=mvsk_objective(lam,mnew)
            if(objnew>obj) then
               if(trim(meth)=='RFPA') then
                  call pgd_backtrack_skew(w,grad,obj,parameters,lam,eta,bet,what,objnew)
               else
                  do i=1,50
                     alpha=0.5_dp*(alpha-1.0_dp)
                     what=w-2.0_dp*alpha*r+alpha*alpha*vvec
                     call evaluate_skew_t_details(what,parameters,mnew,g1)
                     call project_simplex(what-objective_gradient(lam,g1)/tauv,w1)
                     what=w1
                     call evaluate_skew_t_details(what,parameters,mnew)
                     objnew=mvsk_objective(lam,mnew)
                     if(objnew<=obj) exit
                  end do
               end if
            end if
            w=what
         end select
         call project_simplex(w,what)
         w=what
         call evaluate_skew_t_details(w,parameters,mnew,g,h2,h3,h4)
         objnew=mvsk_objective(lam,mnew)
         hist(iter+1)=scale*objnew
         wc=vector_norm2(w-wold)/max(vector_norm2(wold),1.0e-15_dp)<wt
         fc=abs(objnew-obj)<ft
         sc=scale*objnew<=sv
         obj=objnew
         m=mnew
         if(sv<=-0.5_dp*huge(1.0_dp)) then
            if(wc .and. fc) exit
         else
            if(sc) exit
         end if
      end do
      hist(1)=hist(1)*scale
      call finish_result(result,w,m,hist,iter,scale*obj,iter<niter)
   end subroutine design_mvsk_portfolio_via_skew_t

   subroutine design_mvsktilting_portfolio_via_sample_moments(d,statistics,result, &
      w_init,w0,w0_moments,leverage,kappa,method,gamma,zeta,maxiter,ftol,wtol,theta,stopval)
      real(dp), intent(in) :: d(4)
      type(sample_moments), intent(in) :: statistics
      type(portfolio_result), intent(out) :: result
      real(dp), intent(in), optional :: w_init(:),w0(:),w0_moments(4)
      real(dp), intent(in), optional :: leverage,kappa,gamma,zeta,ftol,wtol,theta,stopval
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: w(:),base(:),old(:),cand(:),g(:,:),hist(:),dir(:)
      real(dp) :: m(4),m0(4),imp(4),delta,old_delta,kap,gam,zet,ft,wt,sv,lev,th
      real(dp) :: eta,temp,weights(4),trial_delta,trial_m(4),trial_imp(4),den
      character(len=16) :: meth
      integer :: n,niter,iter,active,j
      logical :: wc,fc,sc

      call clear_portfolio_result(result)
      n=statistics%nassets
      if(n<1 .or. any(d<=0.0_dp)) then
         call fail_result(result,hop_invalid_argument,'invalid tilting inputs')
         return
      end if
      lev=1.0_dp
      if(present(leverage)) lev=leverage
      if(abs(lev-1.0_dp)>1.0e-12_dp) then
         call fail_result(result,hop_invalid_argument,'only leverage=1 is supported')
         return
      end if
      kap=0.0_dp
      if(present(kappa)) kap=max(0.0_dp,kappa)
      meth='Q-MVSKT'
      if(present(method)) meth=to_upper(trim(method))
      if(meth/='Q-MVSKT' .and. meth/='L-MVSKT') then
         call fail_result(result,hop_invalid_argument,'method must be Q-MVSKT or L-MVSKT')
         return
      end if
      gam=1.0_dp
      if(present(gamma)) gam=gamma
      zet=1.0e-8_dp
      if(present(zeta)) zet=zeta
      niter=100
      if(present(maxiter)) niter=maxiter
      ft=1.0e-5_dp
      if(present(ftol)) ft=ftol
      wt=1.0e-5_dp
      if(present(wtol)) wt=wtol
      th=0.5_dp
      if(present(theta)) th=theta
      sv=-huge(1.0_dp)
      if(present(stopval)) sv=stopval

      allocate(w(n),base(n),old(n),cand(n),g(4,n),dir(n),hist(niter+1))
      if(present(w0)) then
         if(size(w0)/=n) then
            call fail_result(result,hop_dimension_mismatch,'w0 has wrong length')
            return
         end if
         call project_simplex(w0,base)
      else if(present(w_init)) then
         if(size(w_init)/=n) then
            call fail_result(result,hop_dimension_mismatch,'w_init has wrong length')
            return
         end if
         call project_simplex(w_init,base)
      else
         base=1.0_dp/real(n,dp)
      end if
      if(present(w_init)) then
         call project_simplex(w_init,w)
      else
         w=base
      end if
      call enforce_tracking(w,base,statistics%covariance,kap)
      if(present(w0_moments)) then
         m0=w0_moments
      else
         call evaluate_sample_details(base,statistics,m0)
      end if
      call evaluate_sample_details(w,statistics,m,g)
      imp=tilting_improvement(m,m0,d)
      delta=minval(imp)
      hist(1)=-delta

      do iter=1,niter
         old=w
         old_delta=delta
         if(meth=='L-MVSKT') then
            active=minloc(imp,dim=1)
            dir=tilting_sign(active)*g(active,:)/d(active)
         else
            temp=max(1.0e-4_dp,0.1_dp/sqrt(real(iter,dp)))
            weights=exp(-(imp-minval(imp))/temp)
            weights=weights/sum(weights)
            dir=0.0_dp
            do j=1,4
               dir=dir+weights(j)*tilting_sign(j)*g(j,:)/d(j)
            end do
         end if
         eta=1.0_dp
         do j=1,60
            call project_simplex(w+eta*dir,cand)
            call enforce_tracking(cand,base,statistics%covariance,kap)
            call evaluate_sample_details(cand,statistics,trial_m)
            trial_imp=tilting_improvement(trial_m,m0,d)
            trial_delta=minval(trial_imp)
            if(trial_delta>=delta-1.0e-12_dp) exit
            eta=eta*th
         end do
         w=w+gam*(cand-w)
         call project_simplex(w,cand)
         w=cand
         call enforce_tracking(w,base,statistics%covariance,kap)
         gam=gam*(1.0_dp-zet*gam)
         call evaluate_sample_details(w,statistics,m,g)
         imp=tilting_improvement(m,m0,d)
         delta=minval(imp)
         hist(iter+1)=-delta
         wc=vector_norm2(w-old)<=wt*max(1.0_dp,vector_norm2(old))
         den=max(1.0_dp,abs(delta))
         fc=abs(delta-old_delta)<=ft*den
         sc=-delta<=sv
         if(wc .or. fc .or. sc) exit
      end do
      call finish_result(result,w,m,hist,iter,-delta,iter<niter)
      result%delta=delta
      result%improvement=imp
   end subroutine design_mvsktilting_portfolio_via_sample_moments

   subroutine pgd_backtrack_skew(w,grad,obj,p,lam,initial_eta,beta,wout,objout)
      real(dp), intent(in) :: w(:),grad(:),obj,lam(4),initial_eta,beta
      type(skew_t_parameters), intent(in) :: p
      real(dp), intent(out) :: wout(:),objout
      real(dp) :: eta,m(4),rhs
      integer :: i
      eta=initial_eta
      do i=1,80
         call project_simplex(w-eta*grad,wout)
         call evaluate_skew_t_details(wout,p,m)
         objout=mvsk_objective(lam,m)
         rhs=obj+dot_product(grad,wout-w)+dot_product(wout-w,wout-w)/(2.0_dp*eta)
         if(objout<=rhs+1.0e-14_dp) return
         eta=eta*beta
      end do
   end subroutine pgd_backtrack_skew

   pure function objective_gradient(lam,g) result(grad)
      real(dp), intent(in) :: lam(4),g(:,:)
      real(dp) :: grad(size(g,2))
      grad=-lam(1)*g(1,:)+lam(2)*g(2,:)-lam(3)*g(3,:)+lam(4)*g(4,:)
   end function objective_gradient

   pure function mvsk_objective(lam,m) result(f)
      real(dp), intent(in) :: lam(4),m(4)
      real(dp) :: f
      f=-lam(1)*m(1)+lam(2)*m(2)-lam(3)*m(3)+lam(4)*m(4)
   end function mvsk_objective

   pure function tilting_improvement(m,m0,d) result(imp)
      real(dp), intent(in) :: m(4),m0(4),d(4)
      real(dp) :: imp(4)
      imp=(m-m0)/d
      imp(2)=-imp(2)
      imp(4)=-imp(4)
   end function tilting_improvement

   pure function tilting_sign(i) result(s)
      integer, intent(in) :: i
      real(dp) :: s
      if(i==1 .or. i==3) then
         s=1.0_dp
      else
         s=-1.0_dp
      end if
   end function tilting_sign

   subroutine enforce_tracking(w,w0,covariance,kappa)
      real(dp), intent(inout) :: w(:)
      real(dp), intent(in) :: w0(:),covariance(:,:),kappa
      real(dp) :: diff(size(w)),q,alpha
      if(kappa<=0.0_dp) then
         w=w0
         return
      end if
      diff=w-w0
      q=dot_product(diff,matmul(covariance,diff))
      if(q>kappa*kappa .and. q>0.0_dp) then
         alpha=kappa/sqrt(q)
         w=w0+alpha*diff
      end if
   end subroutine enforce_tracking

   pure function valid_skew_method(method) result(ok)
      character(len=*), intent(in) :: method
      logical :: ok
      select case(trim(method))
      case('L-MVSK','DC','Q-MVSK','SQUAREM','RFPA','PGD')
         ok=.true.
      case default
         ok=.false.
      end select
   end function valid_skew_method

   pure function to_upper(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i,k
      out=text
      do i=1,len(text)
         k=iachar(out(i:i))
         if(k>=iachar('a') .and. k<=iachar('z')) out(i:i)=achar(k-32)
      end do
   end function to_upper

   subroutine fail_result(result,status,message)
      type(portfolio_result), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status=status
      result%message=message
   end subroutine fail_result

   subroutine finish_result(result,w,m,hist,iter,obj,converged)
      type(portfolio_result), intent(inout) :: result
      real(dp), intent(in) :: w(:),m(4),hist(:),obj
      integer, intent(in) :: iter
      logical, intent(in) :: converged
      integer :: nh
      nh=min(iter+1,size(hist))
      allocate(result%w(size(w)),result%objective_history(nh))
      result%w=w
      result%objective_history=hist(1:nh)
      result%moments=m
      result%iterations=iter
      result%converged=converged
      if(converged) then
         result%status=hop_success
         result%message='success'
      else
         result%status=hop_not_converged
         write(result%message,'(a,es12.4)') 'maximum iterations reached; objective=',obj
      end if
   end subroutine finish_result

end module highorder_optimization
