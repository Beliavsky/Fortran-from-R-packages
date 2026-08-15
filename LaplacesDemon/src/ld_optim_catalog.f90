module ld_optim_catalog
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface, vector_func_iface
use ld_random, only: rand_uniform, rand_normal
use ld_linalg, only: inverse_spd, make_positive_definite, outer_product
use ld_numerics, only: numerical_gradient, numerical_hessian, numerical_jacobian
use ld_optim, only: optim_result_t
implicit none
private
public :: newton_maximize, levenberg_marquardt_maximize, nelder_mead_maximize
public :: conjugate_gradient_maximize, dfp_maximize, lbfgs_maximize
public :: hooke_jeeves_maximize, trust_region_maximize, rprop_maximize
public :: sgd_maximize, spg_maximize, sr1_maximize, pso_maximize
public :: genetic_maximize, soma_maximize, hit_and_run_maximize, bhhh_maximize

contains

subroutine newton_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,info,ls
   real(dp) :: x(size(x0)),g(size(x0)),h(size(x0),size(x0)),a(size(x0),size(x0))
   real(dp) :: step(size(x0)),trial(size(x0)),fx,ft,eps,alpha
   n=size(x0); maxit=200; if(present(max_iter)) maxit=max_iter
   eps=1.0e-7_dp; if(present(tol)) eps=tol; x=x0; fx=f(x)
   do it=1,maxit
      call numerical_gradient(f,x,g); if(maxval(abs(g))<eps) exit
      call numerical_hessian(f,x,h); a=-h; call make_positive_definite(a,info=info)
      if(info/=0) then; step=g; else; call solve_posdef(a,g,step,info); if(info/=0) step=g; end if
      alpha=1.0_dp
      do ls=1,30
         trial=x+alpha*step; ft=f(trial)
         if(ft>=fx+1.0e-4_dp*alpha*dot_product(g,step)) exit
         alpha=0.5_dp*alpha
      end do
      if(alpha<1.0e-12_dp) exit
      x=trial; fx=ft
      if(maxval(abs(alpha*step))<eps*(1.0_dp+maxval(abs(x)))) exit
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine newton_maximize

subroutine levenberg_marquardt_maximize(f,x0,res,max_iter,tol,lambda0)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,lambda0
   integer :: n,it,maxit,info,i
   real(dp) :: x(size(x0)),g(size(x0)),h(size(x0),size(x0)),a(size(x0),size(x0))
   real(dp) :: step(size(x0)),trial(size(x0)),fx,ft,eps,lam
   n=size(x0); maxit=300; if(present(max_iter)) maxit=max_iter
   eps=1.0e-7_dp; if(present(tol)) eps=tol; lam=1.0e-2_dp; if(present(lambda0)) lam=lambda0
   x=x0; fx=f(x)
   do it=1,maxit
      call numerical_gradient(f,x,g); if(maxval(abs(g))<eps) exit
      call numerical_hessian(f,x,h); a=-h
      do i=1,n; a(i,i)=a(i,i)+lam; end do
      call make_positive_definite(a,info=info); call solve_posdef(a,g,step,info)
      if(info/=0) then; lam=lam*10.0_dp; cycle; end if
      trial=x+step; ft=f(trial)
      if(ft>fx) then; x=trial; fx=ft; lam=max(lam/3.0_dp,1.0e-12_dp); else; lam=lam*5.0_dp; end if
      if(maxval(abs(step))<eps*(1.0_dp+maxval(abs(x)))) exit
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine levenberg_marquardt_maximize

subroutine nelder_mead_maximize(f,x0,res,max_iter,tol,initial_step)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,initial_step
   integer :: n,it,maxit,i,lo,hi,second
   real(dp) :: eps,delta,fr,fe,fc
   real(dp), allocatable :: s(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
   n=size(x0); maxit=1000; if(present(max_iter)) maxit=max_iter
   eps=1.0e-8_dp; if(present(tol)) eps=tol; delta=0.1_dp; if(present(initial_step)) delta=initial_step
   allocate(s(n+1,n),fv(n+1),cent(n),xr(n),xe(n),xc(n)); s(1,:)=x0
   do i=1,n; s(i+1,:)=x0; s(i+1,i)=s(i+1,i)+delta*max(1.0_dp,abs(x0(i))); end do
   do i=1,n+1; fv(i)=f(s(i,:)); end do
   do it=1,maxit
      call order_simplex(fv,hi,second,lo)
      if(maxval(abs(fv-fv(lo)))<eps .and. maxval(abs(s-spread(s(lo,:),1,n+1)))<sqrt(eps)) exit
      cent=(sum(s,dim=1)-s(hi,:))/real(n,dp); xr=cent+(cent-s(hi,:)); fr=f(xr)
      if(fr>fv(lo)) then
         xe=cent+2.0_dp*(xr-cent); fe=f(xe)
         if(fe>fr) then; s(hi,:)=xe; fv(hi)=fe; else; s(hi,:)=xr; fv(hi)=fr; end if
      else if(fr>fv(second)) then
         s(hi,:)=xr; fv(hi)=fr
      else
         if(fr>fv(hi)) then; xc=cent+0.5_dp*(xr-cent); else; xc=cent+0.5_dp*(s(hi,:)-cent); end if
         fc=f(xc)
         if(fc>fv(hi)) then
            s(hi,:)=xc; fv(hi)=fc
         else
            do i=1,n+1
               if(i/=lo) then; s(i,:)=s(lo,:)+0.5_dp*(s(i,:)-s(lo,:)); fv(i)=f(s(i,:)); end if
            end do
         end if
      end if
   end do
   lo=maxloc(fv,dim=1); allocate(res%par(n)); res%par=s(lo,:); res%value=fv(lo)
   res%iterations=min(it,maxit); res%converged=(it<=maxit)
end subroutine nelder_mead_maximize

subroutine conjugate_gradient_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,ls
   real(dp) :: x(size(x0)),g(size(x0)),gn(size(x0)),d(size(x0)),trial(size(x0)),fx,ft,beta,alpha,eps,den
   n=size(x0); maxit=500; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   x=x0; fx=f(x); call numerical_gradient(f,x,g); d=g
   do it=1,maxit
      if(maxval(abs(g))<eps) exit; if(dot_product(g,d)<=0.0_dp) d=g
      alpha=1.0_dp
      do ls=1,35
         trial=x+alpha*d; ft=f(trial); if(ft>=fx+1.0e-4_dp*alpha*dot_product(g,d)) exit; alpha=alpha*0.5_dp
      end do
      if(alpha<1.0e-12_dp) exit
      call numerical_gradient(f,trial,gn); den=max(dot_product(g,g),tiny(1.0_dp))
      beta=max(0.0_dp,dot_product(gn,gn-g)/den); d=gn+beta*d; x=trial; g=gn; fx=ft
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine conjugate_gradient_maximize

subroutine dfp_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,ls
   real(dp) :: x(size(x0)),xn(size(x0)),g(size(x0)),gn(size(x0)),d(size(x0))
   real(dp) :: s(size(x0)),y(size(x0)),hy(size(x0)),h(size(x0),size(x0)),fx,fn,alpha,eps,sy,yhy
   n=size(x0); maxit=500; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   h=identity_local(n); x=x0; fx=f(x); call numerical_gradient(f,x,g)
   do it=1,maxit
      if(maxval(abs(g))<eps) exit; d=matmul(h,g); if(dot_product(g,d)<=0.0_dp) d=g
      alpha=1.0_dp
      do ls=1,35; xn=x+alpha*d; fn=f(xn); if(fn>=fx+1.0e-4_dp*alpha*dot_product(g,d)) exit; alpha=alpha*0.5_dp; end do
      if(alpha<1.0e-12_dp) exit
      call numerical_gradient(f,xn,gn); s=xn-x; y=-(gn-g); sy=dot_product(s,y)
      if(sy>1.0e-12_dp) then
         hy=matmul(h,y); yhy=dot_product(y,hy)
         if(yhy>1.0e-12_dp) h=h+outer_product(s,s)/sy-outer_product(hy,hy)/yhy
      end if
      x=xn; g=gn; fx=fn
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine dfp_maximize

subroutine lbfgs_maximize(f,x0,res,max_iter,tol,memory)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter,memory
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,m,used,i,j,ls,idx
   real(dp) :: x(size(x0)),xn(size(x0)),g(size(x0)),gn(size(x0)),q(size(x0)),r(size(x0)),d(size(x0))
   real(dp) :: fx,fn,alpha_ls,eps,gamma0,beta
   real(dp), allocatable :: ss(:,:),yy(:,:),rho(:),aa(:)
   n=size(x0); maxit=500; if(present(max_iter)) maxit=max_iter; m=7; if(present(memory)) m=max(1,memory)
   eps=1.0e-7_dp; if(present(tol)) eps=tol; allocate(ss(m,n),yy(m,n),rho(m),aa(m)); ss=0; yy=0; rho=0
   x=x0; fx=f(x); call numerical_gradient(f,x,g); used=0
   do it=1,maxit
      if(maxval(abs(g))<eps) exit; q=-g
      do j=used,1,-1; idx=1+mod(it-j-1,m); aa(idx)=rho(idx)*dot_product(ss(idx,:),q); q=q-aa(idx)*yy(idx,:); end do
      if(used>0) then
         idx=1+mod(it-2,m)
         gamma0=dot_product(ss(idx,:),yy(idx,:)) / &
              max(dot_product(yy(idx,:),yy(idx,:)),tiny(1.0_dp))
      else
         gamma0=1.0_dp
      end if
      r=gamma0*q
      do j=1,used; idx=1+mod(it-used+j-2,m); beta=rho(idx)*dot_product(yy(idx,:),r); r=r+ss(idx,:)*(aa(idx)-beta); end do
      d=-r; if(dot_product(g,d)<=0.0_dp) d=g
      alpha_ls=1.0_dp
      do ls=1,35; xn=x+alpha_ls*d; fn=f(xn); if(fn>=fx+1.0e-4_dp*alpha_ls*dot_product(g,d)) exit; alpha_ls=0.5_dp*alpha_ls; end do
      if(alpha_ls<1.0e-12_dp) exit; call numerical_gradient(f,xn,gn)
      idx=1+mod(it-1,m); ss(idx,:)=xn-x; yy(idx,:)=-(gn-g)
      if(dot_product(ss(idx,:),yy(idx,:))>1.0e-12_dp) then
         rho(idx)=1.0_dp/dot_product(ss(idx,:),yy(idx,:))
         used=min(m,used+1)
      end if
      x=xn; g=gn; fx=fn
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine lbfgs_maximize

subroutine hooke_jeeves_maximize(f,x0,res,max_iter,tol,step0)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,step0
   integer :: n,it,maxit,j
   real(dp) :: x(size(x0)),base(size(x0)),trial(size(x0)),fx,fb,ft,step,eps
   n=size(x0); maxit=2000; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   step=0.5_dp; if(present(step0)) step=step0; x=x0; fx=f(x)
   do it=1,maxit
      base=x; fb=fx
      do j=1,n
         trial=x; trial(j)=trial(j)+step; ft=f(trial)
         if(ft>fx) then; x=trial; fx=ft; cycle; end if
         trial=x; trial(j)=trial(j)-step; ft=f(trial); if(ft>fx) then; x=trial; fx=ft; end if
      end do
      if(fx>fb) then; trial=x+(x-base); ft=f(trial); if(ft>fx) then; x=trial; fx=ft; end if
      else; step=0.5_dp*step; end if
      if(step<eps*(1.0_dp+maxval(abs(x)))) exit
   end do
   allocate(res%par(n)); res%par=x; res%value=fx; res%iterations=min(it,maxit); res%converged=(step<sqrt(eps))
end subroutine hooke_jeeves_maximize

subroutine trust_region_maximize(f,x0,res,max_iter,tol,radius0)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,radius0
   integer :: n,it,maxit,info
   real(dp) :: x(size(x0)),g(size(x0)),h(size(x0),size(x0)),a(size(x0),size(x0)),step(size(x0)),trial(size(x0))
   real(dp) :: fx,ft,eps,rad,norms,pred,rho
   n=size(x0); maxit=300; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   rad=1.0_dp; if(present(radius0)) rad=radius0; x=x0; fx=f(x)
   do it=1,maxit
      call numerical_gradient(f,x,g); if(maxval(abs(g))<eps) exit; call numerical_hessian(f,x,h); a=-h
      call make_positive_definite(a,info=info); call solve_posdef(a,g,step,info); if(info/=0) step=g
      norms=sqrt(dot_product(step,step)); if(norms>rad) step=step*(rad/norms)
      trial=x+step; ft=f(trial); pred=dot_product(g,step)+0.5_dp*dot_product(step,matmul(h,step))
      if(pred<=0.0_dp) pred=dot_product(g,step); rho=(ft-fx)/max(pred,1.0e-16_dp)
      if(rho>0.1_dp) then; x=trial; fx=ft; end if
      if(rho<0.25_dp) rad=0.25_dp*rad; if(rho>0.75_dp .and. abs(norms-rad)<0.1_dp*rad) rad=min(2.0_dp*rad,100.0_dp)
      if(rad<eps) exit
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine trust_region_maximize

subroutine rprop_maximize(f,x0,res,max_iter,tol,delta0)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,delta0
   integer :: n,it,maxit,j
   real(dp) :: x(size(x0)),g(size(x0)),gold(size(x0)),delta(size(x0)),fx,eps,d0,prod
   n=size(x0); maxit=1000; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   d0=0.1_dp; if(present(delta0)) d0=delta0; x=x0; gold=0.0_dp; delta=d0
   do it=1,maxit
      call numerical_gradient(f,x,g); if(maxval(abs(g))<eps) exit
      do j=1,n
         prod=g(j)*gold(j)
         if(prod>0.0_dp) delta(j)=min(delta(j)*1.2_dp,50.0_dp)
         if(prod<0.0_dp) then; delta(j)=max(delta(j)*0.5_dp,1.0e-8_dp); g(j)=0.0_dp; end if
         if(g(j)>0.0_dp) x(j)=x(j)+delta(j); if(g(j)<0.0_dp) x(j)=x(j)-delta(j)
      end do
      gold=g
   end do
   fx=f(x); call set_result(res,x,fx,it,maxit,eps,f)
end subroutine rprop_maximize

subroutine sgd_maximize(f,x0,res,max_iter,tol,learning_rate,decay)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,learning_rate,decay
   integer :: it,maxit
   real(dp) :: x(size(x0)),g(size(x0)),eps,lr,dec
   maxit=2000; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   lr=0.05_dp; if(present(learning_rate)) lr=learning_rate; dec=0.55_dp; if(present(decay)) dec=decay; x=x0
   do it=1,maxit
      call numerical_gradient(f,x,g); if(maxval(abs(g))<eps) exit; x=x+lr/(real(it,dp)**dec)*g
   end do
   call set_result(res,x,f(x),it,maxit,eps,f)
end subroutine sgd_maximize

subroutine spg_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: it,maxit,ls
   real(dp) :: x(size(x0)),xn(size(x0)),g(size(x0)),gn(size(x0)),s(size(x0)),y(size(x0)),eps,alpha,fx,fn,bb
   maxit=1000; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   x=x0; fx=f(x); call numerical_gradient(f,x,g); bb=1.0_dp
   do it=1,maxit
      if(maxval(abs(g))<eps) exit; alpha=bb
      do ls=1,30; xn=x+alpha*g; fn=f(xn); if(fn>=fx+1.0e-4_dp*alpha*dot_product(g,g)) exit; alpha=0.5_dp*alpha; end do
      if(alpha<1.0e-12_dp) exit; call numerical_gradient(f,xn,gn); s=xn-x; y=gn-g
      if(abs(dot_product(s,y))>1.0e-14_dp) bb=abs(dot_product(s,s)/dot_product(s,y)); bb=min(max(bb,1.0e-6_dp),1.0e3_dp)
      x=xn; g=gn; fx=fn
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine spg_maximize

subroutine sr1_maximize(f,x0,res,max_iter,tol)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,ls
   real(dp) :: x(size(x0)),xn(size(x0)),g(size(x0)),gn(size(x0)),h(size(x0),size(x0)),d(size(x0))
   real(dp) :: s(size(x0)),y(size(x0)),u(size(x0)),eps,alpha,fx,fn,den
   n=size(x0); maxit=500; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   h=identity_local(n); x=x0; fx=f(x); call numerical_gradient(f,x,g)
   do it=1,maxit
      if(maxval(abs(g))<eps) exit; d=matmul(h,g); if(dot_product(g,d)<=0.0_dp) d=g
      alpha=1.0_dp
      do ls=1,35; xn=x+alpha*d; fn=f(xn); if(fn>=fx+1.0e-4_dp*alpha*dot_product(g,d)) exit; alpha=0.5_dp*alpha; end do
      if(alpha<1.0e-12_dp) exit; call numerical_gradient(f,xn,gn); s=xn-x; y=-(gn-g); u=s-matmul(h,y); den=dot_product(u,y)
      if(abs(den)>1.0e-8_dp*sqrt(dot_product(u,u)*dot_product(y,y))) h=h+outer_product(u,u)/den
      x=xn; g=gn; fx=fn
   end do
   call set_result(res,x,fx,it,maxit,eps,f)
end subroutine sr1_maximize

subroutine pso_maximize(f,x0,res,max_iter,tol,n_particles,spread)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter,n_particles
   real(dp), intent(in), optional :: tol,spread
   integer :: n,np,it,maxit,i,j,best,stall
   real(dp) :: eps,spr,w,c1,c2
   real(dp), allocatable :: pos(:,:),vel(:,:),pbest(:,:),pval(:),gpos(:),trial(:)
   n=size(x0); np=max(10,5*n); if(present(n_particles)) np=max(4,n_particles); maxit=500; if(present(max_iter)) maxit=max_iter
   eps=1.0e-6_dp; if(present(tol)) eps=tol; spr=1.0_dp; if(present(spread)) spr=spread; w=0.72_dp; c1=1.49_dp; c2=1.49_dp
   allocate(pos(np,n),vel(np,n),pbest(np,n),pval(np),gpos(n),trial(n))
   do i=1,np
      do j=1,n
         pos(i,j)=x0(j)+spr*rand_normal()
         vel(i,j)=0.1_dp*spr*rand_normal()
      end do
      pbest(i,:)=pos(i,:); pval(i)=f(pos(i,:))
   end do
   best=maxloc(pval,dim=1); gpos=pbest(best,:); stall=0
   do it=1,maxit
      do i=1,np
         do j=1,n
            vel(i,j)=w*vel(i,j)+c1*rand_uniform()*(pbest(i,j)-pos(i,j))+c2*rand_uniform()*(gpos(j)-pos(i,j))
         end do
         pos(i,:)=pos(i,:)+vel(i,:)
         if(f(pos(i,:))>pval(i)) then; pbest(i,:)=pos(i,:); pval(i)=f(pos(i,:)); end if
      end do
      best=maxloc(pval,dim=1); trial=gpos; gpos=pbest(best,:)
      if(maxval(abs(gpos-trial))<eps) then; stall=stall+1; else; stall=0; end if
      if(stall>=30) exit
   end do
   allocate(res%par(n)); res%par=gpos; res%value=f(gpos); res%iterations=min(it,maxit); res%converged=(it<=maxit)
end subroutine pso_maximize

subroutine genetic_maximize(f,x0,res,max_iter,tol,pop_size,spread,mutation_sd)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter,pop_size
   real(dp), intent(in), optional :: tol,spread,mutation_sd
   integer :: n,np,it,maxit,i,j,a,b,best,elite,stall
   real(dp) :: eps,spr,mut,alpha,oldbest,newbest
   real(dp), allocatable :: pop(:,:),next(:,:),val(:)
   n=size(x0); np=max(20,8*n); if(present(pop_size)) np=max(6,pop_size); maxit=400; if(present(max_iter)) maxit=max_iter
   eps=1.0e-6_dp; if(present(tol)) eps=tol
   spr=1.0_dp; if(present(spread)) spr=spread
   mut=0.05_dp*spr; if(present(mutation_sd)) mut=mutation_sd
   allocate(pop(np,n),next(np,n),val(np)); do i=1,np; do j=1,n; pop(i,j)=x0(j)+spr*rand_normal(); end do; val(i)=f(pop(i,:)); end do
   oldbest=-huge(1.0_dp); elite=max(1,np/10); stall=0
   do it=1,maxit
      call sort_population(pop,val)
      next(1:elite,:)=pop(1:elite,:)
      do i=elite+1,np
         a=tournament_index(val); b=tournament_index(val); alpha=rand_uniform(); next(i,:)=alpha*pop(a,:)+(1.0_dp-alpha)*pop(b,:)
         do j=1,n; if(rand_uniform()<0.2_dp) next(i,j)=next(i,j)+mut*rand_normal(); end do
      end do
      pop=next; do i=1,np; val(i)=f(pop(i,:)); end do; newbest=maxval(val)
      if(abs(newbest-oldbest)<eps*(1.0_dp+abs(newbest))) then; stall=stall+1; else; stall=0; end if
      oldbest=newbest; if(stall>=30) exit
   end do
   best=maxloc(val,dim=1); allocate(res%par(n)); res%par=pop(best,:)
   res%value=val(best); res%iterations=min(it,maxit); res%converged=(it<=maxit)
end subroutine genetic_maximize

subroutine soma_maximize(f,x0,res,max_iter,tol,pop_size,spread,path_length,step)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter,pop_size
   real(dp), intent(in), optional :: tol,spread,path_length,step
   integer :: n,np,it,maxit,i,j,best,leader,ns,s
   real(dp) :: eps,spr,path,ds,t,bestv,old
   real(dp), allocatable :: pop(:,:),val(:),cand(:),mask(:)
   n=size(x0); np=max(10,4*n); if(present(pop_size)) np=max(4,pop_size); maxit=300; if(present(max_iter)) maxit=max_iter
   eps=1.0e-6_dp; if(present(tol)) eps=tol
   spr=1.0_dp; if(present(spread)) spr=spread
   path=3.0_dp; if(present(path_length)) path=path_length
   ds=0.3_dp; if(present(step)) ds=step; ns=max(1,ceiling(path/ds)); allocate(pop(np,n),val(np),cand(n),mask(n))
   do i=1,np; do j=1,n; pop(i,j)=x0(j)+spr*rand_normal(); end do; val(i)=f(pop(i,:)); end do; old=-huge(1.0_dp)
   do it=1,maxit
      leader=maxloc(val,dim=1)
      do i=1,np
         if(i==leader) cycle; bestv=val(i); cand=pop(i,:)
         do s=1,ns
            t=min(path,real(s,dp)*ds); do j=1,n; mask(j)=merge(1.0_dp,0.0_dp,rand_uniform()<0.5_dp); end do
            if(sum(mask)==0.0_dp) mask(1+mod(s-1,n))=1.0_dp
            if(f(pop(i,:)+t*mask*(pop(leader,:)-pop(i,:)))>bestv) then
               cand=pop(i,:)+t*mask*(pop(leader,:)-pop(i,:)); bestv=f(cand)
            end if
         end do
         pop(i,:)=cand; val(i)=bestv
      end do
      best=maxloc(val,dim=1); if(abs(val(best)-old)<eps*(1.0_dp+abs(val(best)))) exit; old=val(best)
   end do
   best=maxloc(val,dim=1); allocate(res%par(n)); res%par=pop(best,:)
   res%value=val(best); res%iterations=min(it,maxit); res%converged=(it<=maxit)
end subroutine soma_maximize

subroutine hit_and_run_maximize(f,x0,res,max_iter,tol,step0)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: x0(:)
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol,step0
   integer :: n,it,maxit,j
   real(dp) :: x(size(x0)),dir(size(x0)),trial(size(x0)),fx,ft,eps,step,normd,t
   n=size(x0); maxit=2000; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   step=1.0_dp; if(present(step0)) step=step0; x=x0; fx=f(x)
   do it=1,maxit
      do j=1,n
         dir(j)=rand_normal()
      end do
      normd=sqrt(dot_product(dir,dir)); dir=dir/max(normd,tiny(1.0_dp))
      t=step*rand_normal(); trial=x+t*dir; ft=f(trial)
      if(ft>fx) then; x=trial; fx=ft; else; step=step*0.999_dp; end if; if(step<eps) exit
   end do
   allocate(res%par(n)); res%par=x; res%value=fx; res%iterations=min(it,maxit); res%converged=(step<sqrt(eps))
end subroutine hit_and_run_maximize

subroutine bhhh_maximize(log_components,x0,n_components,res,max_iter,tol)
   procedure(vector_func_iface) :: log_components
   real(dp), intent(in) :: x0(:)
   integer, intent(in) :: n_components
   type(optim_result_t), intent(out) :: res
   integer, intent(in), optional :: max_iter
   real(dp), intent(in), optional :: tol
   integer :: n,it,maxit,info,ls
   real(dp) :: x(size(x0)),trial(size(x0)),step(size(x0)),g(size(x0)),info_mat(size(x0),size(x0)),eps,fx,ft,alpha
   real(dp), allocatable :: c(:),ct(:),jac(:,:)
   n=size(x0); maxit=300; if(present(max_iter)) maxit=max_iter; eps=1.0e-7_dp; if(present(tol)) eps=tol
   allocate(c(n_components),ct(n_components),jac(n_components,n)); x=x0; call log_components(x,c); fx=sum(c)
   do it=1,maxit
      call numerical_jacobian(log_components,x,n_components,jac); g=sum(jac,dim=1); if(maxval(abs(g))<eps) exit
      info_mat=matmul(transpose(jac),jac)+1.0e-10_dp*identity_local(n); call solve_posdef(info_mat,g,step,info); if(info/=0) step=g
      alpha=1.0_dp
      do ls=1,30; trial=x+alpha*step; call log_components(trial,ct); ft=sum(ct); if(ft>=fx) exit; alpha=0.5_dp*alpha; end do
      if(alpha<1.0e-12_dp) exit; x=trial; c=ct; fx=ft
   end do
   allocate(res%par(n)); res%par=x; res%value=fx; res%iterations=min(it,maxit); res%converged=(it<=maxit)
end subroutine bhhh_maximize

subroutine solve_posdef(a,b,x,info)
   use ld_linalg, only: solve_spd
   real(dp), intent(in) :: a(:,:),b(:)
   real(dp), intent(out) :: x(:)
   integer, intent(out) :: info
   call solve_spd(a,b,x,info)
end subroutine solve_posdef

subroutine set_result(res,x,fx,it,maxit,eps,f)
   type(optim_result_t), intent(out) :: res
   real(dp), intent(in) :: x(:),fx,eps
   integer, intent(in) :: it,maxit
   procedure(log_target_iface) :: f
   real(dp) :: g(size(x))
   allocate(res%par(size(x))); res%par=x; res%value=fx; res%iterations=min(it,maxit); call numerical_gradient(f,x,g)
   res%converged=(maxval(abs(g))<sqrt(eps) .or. it<=maxit)
end subroutine set_result

subroutine order_simplex(v,hi,second,lo)
   real(dp), intent(in) :: v(:)
   integer, intent(out) :: hi,second,lo
   integer :: i
   lo=maxloc(v,dim=1); hi=minloc(v,dim=1); second=lo
   do i=1,size(v); if(i/=hi .and. (second==lo .or. v(i)<v(second))) second=i; end do
end subroutine order_simplex

subroutine sort_population(pop,val)
   real(dp), intent(inout) :: pop(:,:),val(:)
   integer :: i,j
   real(dp) :: tv,row(size(pop,2))
   do i=2,size(val)
      tv=val(i); row=pop(i,:); j=i-1
      do while(j>=1)
         if(val(j)>=tv) exit
         val(j+1)=val(j); pop(j+1,:)=pop(j,:); j=j-1
      end do
      val(j+1)=tv; pop(j+1,:)=row
   end do
end subroutine sort_population

integer function tournament_index(val) result(idx)
   real(dp), intent(in) :: val(:)
   integer :: a,b
   a=1+int(rand_uniform()*real(size(val),dp))
   b=1+int(rand_uniform()*real(size(val),dp))
   if(val(a)>=val(b)) then; idx=a; else; idx=b; end if
end function tournament_index

pure function identity_local(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity_local

end module ld_optim_catalog
