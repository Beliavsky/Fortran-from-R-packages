module ld_mcmc_advanced
use ld_kinds, only: dp, pi
use ld_interfaces, only: log_target_iface
use ld_random, only: rand_uniform, rand_normal, rand_student_t, rand_mvn
use ld_linalg, only: chol_lower, make_positive_definite, outer_product, sample_covariance
use ld_numerics, only: numerical_gradient
use ld_mcmc, only: mcmc_result_t
implicit none
private
public :: nuts_sample, hmcda_sample, dram_sample, ram_sample, pcn_sample, twalk_sample

type :: tree_state_t
   real(dp), allocatable :: theta_minus(:), theta_plus(:), r_minus(:), r_plus(:)
   real(dp), allocatable :: grad_minus(:), grad_plus(:), theta_prime(:), grad_prime(:)
   real(dp) :: lp_prime = -huge(1.0_dp)
   real(dp) :: alpha = 0.0_dp
   integer :: nprime = 0, sprime = 0, nalpha = 0
end type tree_state_t

contains

subroutine init_result_adv(res,nkeep,p)
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in) :: nkeep,p
   allocate(res%chain(nkeep,p),res%logp(nkeep))
   res%chain=0.0_dp; res%logp=0.0_dp
   res%acceptance_rate=0.0_dp; res%proposed=0; res%accepted=0
end subroutine init_result_adv

subroutine leapfrog_step(f,theta,r,grad,epsilon,theta_new,r_new,grad_new,lp_new)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: theta(:),r(:),grad(:),epsilon
   real(dp), intent(out) :: theta_new(:),r_new(:),grad_new(:),lp_new
   r_new=r+0.5_dp*epsilon*grad
   theta_new=theta+epsilon*r_new
   lp_new=f(theta_new)
   if (.not. ieee_finite_local(lp_new)) then
      grad_new=0.0_dp
      return
   end if
   call numerical_gradient(f,theta_new,grad_new)
   r_new=r_new+0.5_dp*epsilon*grad_new
end subroutine leapfrog_step

pure logical function ieee_finite_local(x)
   real(dp), intent(in) :: x
   ieee_finite_local=(x==x .and. abs(x)<huge(1.0_dp))
end function ieee_finite_local

pure logical function nuts_stop(theta_minus,theta_plus,r_minus,r_plus)
   real(dp), intent(in) :: theta_minus(:),theta_plus(:),r_minus(:),r_plus(:)
   real(dp) :: d(size(theta_minus))
   d=theta_plus-theta_minus
   nuts_stop=(dot_product(d,r_minus)>=0.0_dp .and. dot_product(d,r_plus)>=0.0_dp)
end function nuts_stop

recursive subroutine build_tree(f,theta,r,grad,logu,v,j,epsilon,joint0,out)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: theta(:),r(:),grad(:),logu,epsilon,joint0
   integer, intent(in) :: v,j
   type(tree_state_t), intent(out) :: out
   type(tree_state_t) :: left,right
   real(dp) :: theta1(size(theta)),r1(size(theta)),grad1(size(theta)),lp1,joint1,den
   integer :: p
   p=size(theta)
   call alloc_tree(out,p)
   if (j==0) then
      call leapfrog_step(f,theta,r,grad,real(v,dp)*epsilon,theta1,r1,grad1,lp1)
      joint1=lp1-0.5_dp*dot_product(r1,r1)
      out%theta_minus=theta1; out%theta_plus=theta1
      out%r_minus=r1; out%r_plus=r1
      out%grad_minus=grad1; out%grad_plus=grad1
      out%theta_prime=theta1; out%grad_prime=grad1; out%lp_prime=lp1
      if (ieee_finite_local(joint1) .and. logu<=joint1) out%nprime=1
      if (ieee_finite_local(joint1) .and. logu<joint1+1000.0_dp) out%sprime=1
      if (ieee_finite_local(joint1)) then
         out%alpha=min(1.0_dp,exp(min(0.0_dp,joint1-joint0)))
      else
         out%alpha=0.0_dp
      end if
      out%nalpha=1
      return
   end if
   call build_tree(f,theta,r,grad,logu,v,j-1,epsilon,joint0,left)
   call copy_tree(left,out)
   if (left%sprime==1) then
      if (v==-1) then
         call build_tree(f,left%theta_minus,left%r_minus,left%grad_minus,logu,v,j-1,epsilon,joint0,right)
         out%theta_minus=right%theta_minus; out%r_minus=right%r_minus; out%grad_minus=right%grad_minus
      else
         call build_tree(f,left%theta_plus,left%r_plus,left%grad_plus,logu,v,j-1,epsilon,joint0,right)
         out%theta_plus=right%theta_plus; out%r_plus=right%r_plus; out%grad_plus=right%grad_plus
      end if
      den=real(max(1,left%nprime+right%nprime),dp)
      if (right%nprime>0 .and. rand_uniform()<real(right%nprime,dp)/den) then
         out%theta_prime=right%theta_prime; out%grad_prime=right%grad_prime
         out%lp_prime=right%lp_prime
      end if
      out%nprime=left%nprime+right%nprime
      if (right%sprime==1 .and. nuts_stop(out%theta_minus,out%theta_plus,out%r_minus,out%r_plus)) then
         out%sprime=1
      else
         out%sprime=0
      end if
      out%alpha=left%alpha+right%alpha
      out%nalpha=left%nalpha+right%nalpha
   end if
end subroutine build_tree

subroutine alloc_tree(t,p)
   type(tree_state_t), intent(inout) :: t
   integer, intent(in) :: p
   if (.not.allocated(t%theta_minus)) then
      allocate(t%theta_minus(p),t%theta_plus(p),t%r_minus(p),t%r_plus(p))
      allocate(t%grad_minus(p),t%grad_plus(p),t%theta_prime(p),t%grad_prime(p))
   end if
end subroutine alloc_tree

subroutine copy_tree(a,b)
   type(tree_state_t), intent(in) :: a
   type(tree_state_t), intent(inout) :: b
   call alloc_tree(b,size(a%theta_minus))
   b%theta_minus=a%theta_minus; b%theta_plus=a%theta_plus
   b%r_minus=a%r_minus; b%r_plus=a%r_plus
   b%grad_minus=a%grad_minus; b%grad_plus=a%grad_plus
   b%theta_prime=a%theta_prime; b%grad_prime=a%grad_prime
   b%lp_prime=a%lp_prime; b%alpha=a%alpha
   b%nprime=a%nprime; b%sprime=a%sprime; b%nalpha=a%nalpha
end subroutine copy_tree

subroutine reasonable_epsilon(f,x,grad,epsilon)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x(:),grad(:)
   real(dp), intent(out) :: epsilon
   real(dp) :: r(size(x)),xn(size(x)),rn(size(x)),gn(size(x)),lp,lpn,loga
   integer :: i,a,k
   do i=1,size(x); r(i)=rand_normal(); end do
   lp=f(x); epsilon=1.0_dp
   call leapfrog_step(f,x,r,grad,epsilon,xn,rn,gn,lpn)
   loga=lpn-lp-0.5_dp*(dot_product(rn,rn)-dot_product(r,r))
   if (.not.ieee_finite_local(loga)) loga=-huge(1.0_dp)
   if (loga>log(0.5_dp)) then; a=1; else; a=-1; end if
   do k=1,30
      if ((a==1 .and. loga<=log(0.5_dp)) .or. (a==-1 .and. loga>=log(0.5_dp))) exit
      if (a==1) then; epsilon=epsilon*2.0_dp; else; epsilon=epsilon*0.5_dp; end if
      call leapfrog_step(f,x,r,grad,epsilon,xn,rn,gn,lpn)
      loga=lpn-lp-0.5_dp*(dot_product(rn,rn)-dot_product(r,r))
      if (.not.ieee_finite_local(loga)) loga=-huge(1.0_dp)
   end do
   epsilon=max(epsilon,1.0e-6_dp)
end subroutine reasonable_epsilon

subroutine nuts_sample(f,x0,n_iter,burn,thin,res,adapt_steps,target_accept,max_depth,initial_step,final_step)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: adapt_steps,max_depth
   real(dp), intent(in), optional :: target_accept,initial_step
   real(dp), intent(out), optional :: final_step
   integer :: p,it,k,nkeep,a_steps,depth,v,n,nold,s,maxd,i
   real(dp) :: x(size(x0)),grad(size(x0)),r0(size(x0)),lp,joint,logu,eps,delta
   real(dp) :: epsbar,hbar,mu,gamma,t0,kappa,eta,accstat
   real(dp) :: thminus(size(x0)),thplus(size(x0)),rminus(size(x0)),rplus(size(x0))
   real(dp) :: grminus(size(x0)),grplus(size(x0))
   type(tree_state_t) :: tree
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   a_steps=burn; if(present(adapt_steps)) a_steps=max(0,adapt_steps)
   maxd=10; if(present(max_depth)) maxd=max(1,max_depth)
   delta=0.65_dp; if(present(target_accept)) delta=target_accept
   x=x0; lp=f(x); call numerical_gradient(f,x,grad)
   if(present(initial_step)) then; eps=initial_step; else; call reasonable_epsilon(f,x,grad,eps); end if
   epsbar=1.0_dp; hbar=0.0_dp; gamma=0.05_dp; t0=10.0_dp; kappa=0.75_dp; mu=log(10.0_dp*eps)
   k=0
   do it=1,n_iter
      do i=1,p; r0(i)=rand_normal(); end do
      joint=lp-0.5_dp*dot_product(r0,r0); logu=joint+log(rand_uniform())
      thminus=x; thplus=x; rminus=r0; rplus=r0; grminus=grad; grplus=grad
      n=1; s=1; depth=0; accstat=0.0_dp
      do while(s==1 .and. depth<maxd)
         if(rand_uniform()<0.5_dp) then; v=-1; else; v=1; end if
         if(v==-1) then
            call build_tree(f,thminus,rminus,grminus,logu,v,depth,eps,joint,tree)
            thminus=tree%theta_minus; rminus=tree%r_minus; grminus=tree%grad_minus
         else
            call build_tree(f,thplus,rplus,grplus,logu,v,depth,eps,joint,tree)
            thplus=tree%theta_plus; rplus=tree%r_plus; grplus=tree%grad_plus
         end if
         nold=n
         if(tree%sprime==1 .and. tree%nprime>0) then
            if(rand_uniform()<real(tree%nprime,dp)/real(max(1,nold+tree%nprime),dp)) then
               x=tree%theta_prime; grad=tree%grad_prime; lp=tree%lp_prime
               res%accepted=res%accepted+1
            end if
         end if
         n=n+tree%nprime
         if(tree%sprime==1 .and. nuts_stop(thminus,thplus,rminus,rplus)) then; s=1; else; s=0; end if
         depth=depth+1; res%proposed=res%proposed+1
         if(tree%nalpha>0) accstat=tree%alpha/real(tree%nalpha,dp)
      end do
      if(it<=a_steps .and. a_steps>0) then
         eta=1.0_dp/(real(it,dp)+t0)
         hbar=(1.0_dp-eta)*hbar+eta*(delta-accstat)
         eps=exp(mu-sqrt(real(it,dp))/gamma*hbar)
         eta=real(it,dp)**(-kappa)
         epsbar=exp(eta*log(eps)+(1.0_dp-eta)*log(epsbar))
      else if(it==a_steps+1 .and. a_steps>0) then
         eps=epsbar
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_step)) final_step=eps
end subroutine nuts_sample

subroutine hmcda_sample(f,x0,n_iter,burn,thin,path_length,res,adapt_steps,target_accept,initial_step,final_step,max_leapfrog)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),path_length
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: adapt_steps,max_leapfrog
   real(dp), intent(in), optional :: target_accept,initial_step
   real(dp), intent(out), optional :: final_step
   integer :: p,it,k,nkeep,a_steps,l,nlf,maxlf,i
   real(dp) :: x(size(x0)),q(size(x0)),r(size(x0)),r0(size(x0)),g(size(x0)),lp,lpp
   real(dp) :: eps,delta,h0,h1,alpha,epsbar,hbar,mu,gamma,t0,kappa,eta
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   a_steps=burn; if(present(adapt_steps)) a_steps=max(0,adapt_steps)
   maxlf=1000; if(present(max_leapfrog)) maxlf=max(1,max_leapfrog)
   delta=0.65_dp; if(present(target_accept)) delta=target_accept
   x=x0; lp=f(x); call numerical_gradient(f,x,g)
   if(present(initial_step)) then; eps=initial_step; else; call reasonable_epsilon(f,x,g,eps); end if
   epsbar=1.0_dp; hbar=0.0_dp; gamma=0.05_dp; t0=10.0_dp; kappa=0.75_dp; mu=log(10.0_dp*eps); k=0
   do it=1,n_iter
      do i=1,p; r0(i)=rand_normal(); end do
      q=x; r=r0; h0=-lp+0.5_dp*dot_product(r0,r0); lpp=lp
      nlf=min(maxlf,max(1,nint(path_length/max(eps,1.0e-12_dp))))
      call numerical_gradient(f,q,g); r=r+0.5_dp*eps*g
      do l=1,nlf
         q=q+eps*r; lpp=f(q)
         if(.not.ieee_finite_local(lpp)) exit
         call numerical_gradient(f,q,g)
         if(l<nlf) r=r+eps*g
      end do
      if(ieee_finite_local(lpp)) r=r+0.5_dp*eps*g
      h1=-lpp+0.5_dp*dot_product(r,r)
      if(ieee_finite_local(h1)) then; alpha=min(1.0_dp,exp(min(0.0_dp,h0-h1))); else; alpha=0.0_dp; end if
      res%proposed=res%proposed+1
      if(rand_uniform()<alpha) then
         x=q; lp=lpp; res%accepted=res%accepted+1
      end if
      if(it<=a_steps .and. a_steps>0) then
         eta=1.0_dp/(real(it,dp)+t0); hbar=(1.0_dp-eta)*hbar+eta*(delta-alpha)
         eps=exp(mu-sqrt(real(it,dp))/gamma*hbar); eta=real(it,dp)**(-kappa)
         epsbar=exp(eta*log(eps)+(1.0_dp-eta)*log(epsbar))
      else if(it==a_steps+1 .and. a_steps>0) then
         eps=epsbar
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_step)) final_step=eps
end subroutine hmcda_sample

subroutine dram_sample(f,x0,initial_cov,n_iter,burn,thin,res,adapt_start,periodicity,scale_factor,final_cov)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),initial_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer, intent(in), optional :: adapt_start,periodicity
   real(dp), intent(in), optional :: scale_factor
   real(dp), allocatable, intent(out), optional :: final_cov(:,:)
   integer :: p,it,k,nkeep,astart,period,info
   real(dp) :: x(size(x0)),y1(size(x0)),y2(size(x0)),step(size(x0)),lp,lp1,lp2
   real(dp) :: loga1,loga2,sc,alpha_y2_y1,alpha_x_y1
   real(dp), allocatable :: cov(:,:),cov2(:,:),history(:,:)
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   astart=max(20,2*p); if(present(adapt_start)) astart=adapt_start
   period=20; if(present(periodicity)) period=max(1,periodicity)
   sc=2.38_dp**2/real(p,dp); if(present(scale_factor)) sc=scale_factor
   allocate(cov(p,p),cov2(p,p),history(n_iter,p)); cov=initial_cov; call make_positive_definite(cov); cov2=0.5_dp*cov
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      call rand_mvn(0.0_dp*x,cov,step,info); y1=x+step; lp1=f(y1); res%proposed=res%proposed+1
      loga1=min(0.0_dp,lp1-lp)
      if(log(rand_uniform())<loga1) then
         x=y1; lp=lp1; res%accepted=res%accepted+1
      else
         call rand_mvn(0.0_dp*x,cov2,step,info); y2=x+step; lp2=f(y2); res%proposed=res%proposed+1
         alpha_y2_y1=min(1.0_dp,exp(min(0.0_dp,lp1-lp2)))
         alpha_x_y1=min(1.0_dp,exp(min(0.0_dp,lp1-lp)))
         if(alpha_y2_y1>=1.0_dp) then
            loga2=-huge(1.0_dp)
         else if(alpha_x_y1>=1.0_dp) then
            loga2=0.0_dp
         else
            loga2=min(0.0_dp,lp2-lp+log(max(1.0e-300_dp,1.0_dp-alpha_y2_y1)) &
                 -log(max(1.0e-300_dp,1.0_dp-alpha_x_y1)))
         end if
         if(log(rand_uniform())<loga2) then
            x=y2; lp=lp2; res%accepted=res%accepted+1
         end if
      end if
      history(it,:)=x
      if(it>=astart .and. mod(it,period)==0 .and. it>2) then
         call sample_covariance(history(1:it,:),cov); cov=sc*cov
         cov=cov+identity_local(p)*(sc*1.0e-5_dp); call make_positive_definite(cov); cov2=0.5_dp*cov
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_cov)) then; allocate(final_cov(p,p)); final_cov=cov; end if
end subroutine dram_sample

subroutine ram_sample(f,x0,initial_cov,n_iter,burn,thin,res,target_accept,gamma,n0,student_t_proposal,final_cov)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),initial_cov(:,:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: target_accept,gamma,n0
   logical, intent(in), optional :: student_t_proposal
   real(dp), allocatable, intent(out), optional :: final_cov(:,:)
   integer :: p,it,k,nkeep,i,info
   real(dp) :: x(size(x0)),prop(size(x0)),u(size(x0)),lp,lpp,loga,astar,gam,ninit,eta,den
   real(dp), allocatable :: cov(:,:),l(:,:),a(:,:),covtest(:,:)
   logical :: use_t
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   astar=0.234_dp; if(present(target_accept)) astar=target_accept
   gam=0.66_dp; if(present(gamma)) gam=gamma
   ninit=1.0_dp; if(present(n0)) ninit=n0
   use_t=.false.; if(present(student_t_proposal)) use_t=student_t_proposal
   allocate(cov(p,p),l(p,p),a(p,p),covtest(p,p)); cov=initial_cov; call make_positive_definite(cov); call chol_lower(cov,l,info)
   x=x0; lp=f(x); k=0
   do it=1,n_iter
      do i=1,p
         if(use_t) then; u(i)=rand_student_t(5.0_dp); else; u(i)=rand_normal(); end if
      end do
      prop=x+matmul(l,u); lpp=f(prop); loga=min(0.0_dp,lpp-lp); res%proposed=res%proposed+1
      if(log(rand_uniform())<loga) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      eta=min(1.0_dp,real(p,dp)*(ninit+real(it,dp))**(-gam)); den=max(dot_product(u,u),1.0e-14_dp)
      a=identity_local(p)+eta*(min(1.0_dp,exp(loga))-astar)*outer_product(u,u)/den
      covtest=matmul(l,matmul(a,transpose(l))); call make_positive_definite(covtest)
      call chol_lower(covtest,l,info); if(info==0) cov=covtest
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
   if(present(final_cov)) then; allocate(final_cov(p,p)); final_cov=cov; end if
end subroutine ram_sample

subroutine pcn_sample(f,x0,proposal_cov,beta,n_iter,burn,thin,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),proposal_cov(:,:),beta
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   integer :: p,it,k,nkeep,info
   real(dp) :: x(size(x0)),prop(size(x0)),z(size(x0)),lp,lpp,b
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   b=min(max(beta,1.0e-8_dp),1.0_dp); x=x0; lp=f(x); k=0
   do it=1,n_iter
      call rand_mvn(0.0_dp*x,proposal_cov,z,info)
      prop=sqrt(max(0.0_dp,1.0_dp-b*b))*x+b*z; lpp=f(prop); res%proposed=res%proposed+1
      if(log(rand_uniform())<lpp-lp) then; x=prop; lp=lpp; res%accepted=res%accepted+1; end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lp
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine pcn_sample

subroutine twalk_sample(f,x0,xp0,n_iter,burn,thin,res,at,aw,pphi)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:),xp0(:)
   integer, intent(in) :: n_iter,burn,thin
   type(mcmc_result_t), intent(out) :: res
   real(dp), intent(in), optional :: at,aw,pphi
   integer :: p,it,k,nkeep,i,nphi
   real(dp) :: x(size(x0)),xp(size(x0)),prop(size(x0)),pivot(size(x0)),lpx,lpxp,lpp
   real(dp) :: atv,awv,pp,ker,beta,z,sigma,loga,u,w1,w2
   logical :: move_x,mask(size(x0))
   p=size(x0); nkeep=max(0,(n_iter-burn)/max(1,thin)); call init_result_adv(res,nkeep,p)
   atv=6.0_dp; if(present(at)) atv=at
   awv=1.5_dp; if(present(aw)) awv=aw
   pp=min(1.0_dp,max(1.0_dp/real(max(1,p),dp),4.0_dp/real(max(1,p),dp)))
   if(present(pphi)) pp=min(1.0_dp,max(0.0_dp,pphi))
   x=x0; xp=xp0; lpx=f(x); lpxp=f(xp); k=0
   do it=1,n_iter
      move_x=(rand_uniform()<0.5_dp)
      if(move_x) then; prop=x; pivot=xp; else; prop=xp; pivot=x; end if
      do i=1,p; mask(i)=rand_uniform()<pp; end do
      if(.not.any(mask)) mask(1+int(rand_uniform()*real(p,dp)))=.true.
      nphi=count(mask); ker=rand_uniform(); loga=-huge(1.0_dp)
      if(ker<0.4918_dp) then
         if(rand_uniform()<(atv-1.0_dp)/(2.0_dp*atv)) then
            beta=rand_uniform()**(1.0_dp/(atv+1.0_dp))
         else
            beta=rand_uniform()**(1.0_dp/(1.0_dp-atv))
         end if
         do i=1,p
            if(mask(i)) prop(i)=pivot(i)+beta*(pivot(i)-prop(i))
         end do
         lpp=f(prop)
         if(move_x) then; loga=lpp-lpx+real(nphi-2,dp)*log(beta); else; loga=lpp-lpxp+real(nphi-2,dp)*log(beta); end if
      else if(ker<0.9836_dp) then
         do i=1,p
            if(mask(i)) then
               u=rand_uniform(); z=(awv/(1.0_dp+awv))*(awv*u*u+2.0_dp*u-1.0_dp)
               prop(i)=prop(i)+(prop(i)-pivot(i))*z
            end if
         end do
         lpp=f(prop); if(move_x) then; loga=lpp-lpx; else; loga=lpp-lpxp; end if
      else if(ker<0.9918_dp) then
         sigma=0.0_dp
         do i=1,p; if(mask(i)) sigma=max(sigma,abs(pivot(i)-prop(i))); end do
         sigma=max(sigma,1.0e-12_dp)
         do i=1,p; if(mask(i)) prop(i)=prop(i)+sigma*rand_normal(); end do
         lpp=f(prop)
         w1=twalk_g3(nphi,sigma,prop,pivot); w2=twalk_g3(nphi,sigma,merge(x,xp,move_x),pivot)
         if(move_x) then; loga=lpp-lpx+w1-w2; else; loga=lpp-lpxp+w1-w2; end if
      else
         sigma=0.0_dp
         do i=1,p; if(mask(i)) sigma=max(sigma,abs(pivot(i)-prop(i)))/3.0_dp; end do
         sigma=max(sigma,1.0e-12_dp)
         do i=1,p; if(mask(i)) prop(i)=pivot(i)+sigma*rand_normal(); end do
         lpp=f(prop); if(move_x) then; loga=lpp-lpx; else; loga=lpp-lpxp; end if
      end if
      res%proposed=res%proposed+1
      if(ieee_finite_local(lpp) .and. log(rand_uniform())<min(0.0_dp,loga)) then
         if(move_x) then; x=prop; lpx=lpp; else; xp=prop; lpxp=lpp; end if
         res%accepted=res%accepted+1
      end if
      if(it>burn .and. mod(it-burn,max(1,thin))==0 .and. k<nkeep) then
         k=k+1; res%chain(k,:)=x; res%logp(k)=lpx
      end if
   end do
   if(res%proposed>0) res%acceptance_rate=real(res%accepted,dp)/real(res%proposed,dp)
end subroutine twalk_sample

pure function twalk_g3(nphi,sigma,h,pivot) result(v)
   integer, intent(in) :: nphi
   real(dp), intent(in) :: sigma,h(:),pivot(:)
   real(dp) :: v
   if(nphi<=0) then; v=0.0_dp; else
      v=0.5_dp*real(nphi,dp)*log(2.0_dp*pi)+real(nphi,dp)*log(sigma)+0.5_dp*sum((h-pivot)**2)/(sigma*sigma)
   end if
end function twalk_g3

pure function identity_local(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp
   do i=1,n; a(i,i)=1.0_dp; end do
end function identity_local

end module ld_mcmc_advanced
