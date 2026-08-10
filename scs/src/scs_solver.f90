! SPDX-License-Identifier: GPL-3.0-only
! Homogeneous self-dual Douglas-Rachford SCS iteration translated from
! upstream SCS 3.x (MIT license) and exposed in the R package scs (GPL-3).
module scs_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, ieee_is_nan
   use scs_kinds, only : dp, i4
   use scs_types
   use scs_sparse, only : validate_matrix, accum_by_a, accum_by_atrans, accum_by_p
   use scs_linalg, only : norm_2, norm_inf, safe_div_pos
   use scs_cones, only : validate_cone, set_r_y, project_dual_cone, scale_box_cone
   use scs_normalize, only : normalize_a_p, normalize_b_c
   use scs_ldlt, only : scs_ldlt_factor
   use scs_acceleration, only : aa_workspace
   implicit none
   private
   public :: scs, scs_set_default_settings

   real(dp), parameter :: infeas_negativity_tol=1.0e-9_dp
   integer, parameter :: feasible_iters=1, rescaling_min_iters=100, converged_interval=25
   real(dp), parameter :: tau_factor=10.0_dp, iterate_norm=1.0_dp
   real(dp), parameter :: max_scale_value=1.0e6_dp, min_scale_value=1.0e-6_dp

   type :: residual_state
      integer :: last_iter=-1
      real(dp) :: xt_p_x=0.0_dp,xt_p_x_tau=0.0_dp,ctx=0.0_dp,ctx_tau=0.0_dp,bty=0.0_dp,bty_tau=0.0_dp
      real(dp) :: pobj=0.0_dp,dobj=0.0_dp,gap=0.0_dp,tau=0.0_dp,kap=0.0_dp
      real(dp) :: res_pri=0.0_dp,res_dual=0.0_dp,res_infeas=0.0_dp,res_unbdd_p=0.0_dp,res_unbdd_a=0.0_dp
      real(dp),allocatable :: ax(:),ax_s(:),px(:),aty(:),ax_s_btau(:),px_aty_ctau(:)
      real(dp),allocatable :: x(:),y(:),s(:)
   end type residual_state

contains

   subroutine scs_set_default_settings(stgs)
      type(scs_settings),intent(out)::stgs
      stgs=scs_settings()
   end subroutine scs_set_default_settings

   subroutine scs(data,cone,stgs,sol,info)
      type(scs_data),intent(in)::data
      type(scs_cone),intent(in)::cone
      type(scs_settings),intent(in),optional::stgs
      type(scs_solution),intent(inout)::sol
      type(scs_info),intent(out)::info
      type(scs_settings)::sset
      type(scs_data)::d,dorig
      type(scs_cone)::k
      type(scs_scaling)::scal
      type(scs_ldlt_factor)::fac
      type(aa_workspace)::aa
      type(residual_state)::res
      real(dp),allocatable::u(:),ut(:),v(:),vprev(:),rsk(:),h(:),g(:),diag_r(:)
      real(dp)::box_t,nmb,nmc,sum_log_scale_factor,aa_norm
      real(dp)::t0,t1,solve_start,lin_total,cone_total
      integer::n,m,l,iter,last_scale_update,n_log_scale_factor,scale_updates
      integer::rejected_accel_steps,accepted_accel_steps,aa_status
      logical::ok,have_scal,time_limit_reached

      sset=scs_settings();if(present(stgs))sset=stgs
      call cpu_time(t0)
      info=scs_info()
      if(.not.validate_inputs(data,cone,sset))then
         call failure_solution(data,sol,info);return
      end if
      dorig=data;d=data;k=cone;n=d%n;m=d%m;l=n+m+1
      nmb=norm_inf(dorig%b);nmc=norm_inf(dorig%c)
      have_scal=sset%normalize
      if(have_scal)then
         call normalize_a_p(d%P,d%has_p,d%A,k,scal)
         call normalize_b_c(scal,d%b,d%c)
         call scale_box_cone(k,scal)
      else
         call scale_box_cone(k)
      end if
      allocate(u(l),ut(l),v(l),vprev(l),rsk(l),h(l-1),g(l-1),diag_r(l))
      allocate(res%ax(m),res%ax_s(m),res%ax_s_btau(m),res%px(n),res%aty(n),res%px_aty_ctau(n))
      allocate(res%x(n),res%y(m),res%s(m))
      call set_diag_r(k,sset,diag_r,n,m)
      call fac%factorize(d%A,d%P,d%has_p,diag_r,ok)
      if(.not.ok)then;call failure_solution(data,sol,info);info%status='failure: KKT factorization';return;end if
      h(1:n)=d%c;h(n+1:n+m)=d%b
      call update_g(fac,h,g,n,m)
      if(sset%warm_start .and. allocated(sol%x) .and. allocated(sol%y) .and. allocated(sol%s))then
         call warm_start(sol,scal,have_scal,diag_r,v,n,m)
      else
         v=0.0_dp;v(l)=1.0_dp
      end if
      box_t=1.0_dp;sum_log_scale_factor=0.0_dp;n_log_scale_factor=0;last_scale_update=0;scale_updates=0
      rejected_accel_steps=0;accepted_accel_steps=0
      if(sset%acceleration_lookback/=0) call aa%init(l,abs(sset%acceleration_lookback),sset%acceleration_lookback>0)
      lin_total=0.0_dp;cone_total=0.0_dp;aa_norm=0.0_dp;time_limit_reached=.false.
      call cpu_time(t1);info%setup_time=1000.0_dp*(t1-t0);call cpu_time(solve_start)

      do iter=0,sset%max_iters-1
         aa_norm=0.0_dp
         if(sset%acceleration_lookback/=0 .and. iter>0 .and. mod(iter,sset%acceleration_interval)==0) then
            aa_norm=aa%apply(v,vprev)
         end if
         if(iter>=feasible_iters)call normalize_iterate(v)
         vprev=v
         call cpu_time(t0)
         call project_lin_sys(fac,g,diag_r,v,ut,n,m,iter)
         call cpu_time(t1);lin_total=lin_total+1000.0_dp*(t1-t0)
         call cpu_time(t0)
         u=2.0_dp*ut-v
         call project_dual_cone(u(n+1:n+m),k,diag_r(n+1:n+m),box_t)
         if(iter<feasible_iters)then;u(l)=1.0_dp;else;u(l)=max(u(l),0.0_dp);end if
         call cpu_time(t1);cone_total=cone_total+1000.0_dp*(t1-t0)
         rsk=(v+u-2.0_dp*ut)*diag_r

         if(mod(iter,converged_interval)==0)then
            call populate_residual(data,scal,have_scal,u,rsk,res,iter,n,m)
            info%status_val=convergence_status(res,nmb,nmc,sset)
            if(info%status_val/=scs_unfinished)exit
            if(sset%time_limit_secs>0.0_dp)then
               call cpu_time(t1)
               if(t1-solve_start>sset%time_limit_secs)then;time_limit_reached=.true.;exit;end if
            end if
         end if

         if(sset%verbose .and. mod(iter,250)==0)then
            if(res%last_iter/=iter)call populate_residual(data,scal,have_scal,u,rsk,res,iter,n,m)
            write(*,'(i7,4(1x,es10.3),1x,es10.3)')iter,res%res_pri,res%res_dual,res%gap,0.5_dp*(res%pobj+res%dobj),sset%scale
         end if

         if(sset%adaptive_scale .and. res%last_iter==iter)then
            call maybe_update_scale(d,k,sset,fac,diag_r,g,h,u,ut,v,rsk,res,nmb,nmc, &
               iter,last_scale_update,sum_log_scale_factor,n_log_scale_factor,scale_updates,aa)
         end if
         v=v+sset%alpha*(u-ut)
         if(sset%acceleration_lookback/=0 .and. mod(iter,sset%acceleration_interval)==0 .and. aa_norm>0.0_dp) then
            aa_status=aa%safeguard(v,vprev)
            if(aa_status<0)then;rejected_accel_steps=rejected_accel_steps+1;else;accepted_accel_steps=accepted_accel_steps+1;end if
         end if
      end do

      if(res%last_iter/=iter)call populate_residual(data,scal,have_scal,u,rsk,res,iter,n,m)
      call finalize_solution(res,info%status_val,time_limit_reached,sol,info)
      info%iter=int(iter,i4);info%scale=sset%scale;info%scale_updates=int(scale_updates,i4)
      info%lin_sys_time=lin_total;info%cone_time=cone_total;info%accel_time=0.0_dp
      info%rejected_accel_steps=int(rejected_accel_steps,i4);info%accepted_accel_steps=int(accepted_accel_steps,i4)
      info%kkt_nnz=fac%kkt_nnz;info%factor_nnz=fac%factor_nnz
      info%factorizations=fac%factorizations;info%symbolic_analyses=fac%symbolic_analyses
      info%comp_slack=abs(dot_product(res%s,res%y))
      call cpu_time(t1);info%solve_time=1000.0_dp*(t1-solve_start)
      if(sset%verbose)then
         write(*,'(a,a)')'status: ',trim(info%status)
         write(*,'(a,es14.6)')'objective: ',0.5_dp*(info%pobj+info%dobj)
      end if
   end subroutine scs

   logical function validate_inputs(d,k,s) result(ok)
      type(scs_data),intent(in)::d;type(scs_cone),intent(in)::k;type(scs_settings),intent(in)::s
      ok=.false.
      if(d%m<=0 .or. d%n<=0)return
      if(size(d%b)/=d%m .or. size(d%c)/=d%n)return
      if(d%A%m/=d%m .or. d%A%n/=d%n .or. .not.validate_matrix(d%A))return
      if(d%has_p)then;if(d%P%n/=d%n .or. .not.validate_matrix(d%P,.true.))return;end if
      if(.not.validate_cone(k,d%m))return
      if(s%max_iters<=0 .or. s%eps_abs<0 .or. s%eps_rel<0 .or. s%eps_infeas<0)return
      if(s%alpha<=0 .or. s%alpha>=2 .or. s%rho_x<=0 .or. s%scale<=0)return
      if(s%acceleration_interval<=0)return
      ok=.true.
   end function validate_inputs

   subroutine failure_solution(d,sol,info)
      type(scs_data),intent(in)::d;type(scs_solution),intent(inout)::sol;type(scs_info),intent(inout)::info
      real(dp)::nan
      nan=ieee_value(1.0_dp,ieee_quiet_nan)
      if(allocated(sol%x))deallocate(sol%x);if(allocated(sol%y))deallocate(sol%y);if(allocated(sol%s))deallocate(sol%s)
      allocate(sol%x(max(0,d%n)),sol%y(max(0,d%m)),sol%s(max(0,d%m)));sol%x=nan;sol%y=nan;sol%s=nan
      info%status='failure';info%status_val=scs_failed;info%pobj=nan;info%dobj=nan;info%gap=nan;info%res_pri=nan;info%res_dual=nan
   end subroutine failure_solution

   subroutine set_diag_r(k,s,diag,n,m)
      type(scs_cone),intent(in)::k;type(scs_settings),intent(in)::s;real(dp),intent(out)::diag(:);integer,intent(in)::n,m
      diag(1:n)=s%rho_x;call set_r_y(k,s%scale,diag(n+1:n+m));diag(n+m+1)=tau_factor
   end subroutine set_diag_r

   subroutine update_g(fac,h,g,n,m)
      type(scs_ldlt_factor),intent(in)::fac;real(dp),intent(in)::h(:);real(dp),intent(out)::g(:);integer,intent(in)::n,m
      g=h;g(n+1:n+m)=-g(n+1:n+m);call fac%solve(g)
   end subroutine update_g

   subroutine warm_start(sol,scal,have_scal,diag,v,n,m)
      type(scs_solution),intent(in)::sol;type(scs_scaling),intent(in)::scal;logical,intent(in)::have_scal
      real(dp),intent(in)::diag(:);real(dp),intent(out)::v(:);integer,intent(in)::n,m
      real(dp),allocatable::x(:),y(:),ss(:)
      allocate(x(n),y(m),ss(m));x=sol%x;y=sol%y;ss=sol%s
      if(have_scal)then;x=x/(scal%E/scal%dual_scale);y=y/(scal%D/scal%primal_scale);ss=ss*(scal%D*scal%dual_scale);end if
      where(ieee_is_nan(x))x=0.0_dp;where(ieee_is_nan(y))y=0.0_dp;where(ieee_is_nan(ss))ss=0.0_dp
      v(1:n)=x;v(n+1:n+m)=y+ss/diag(n+1:n+m);v(n+m+1)=1.0_dp
   end subroutine warm_start

   pure real(dp) function dot_r(diag,x,y,nm) result(v)
      real(dp),intent(in)::diag(:),x(:),y(:);integer,intent(in)::nm
      v=sum(x(1:nm)*y(1:nm)*diag(1:nm))
   end function dot_r

   real(dp) function root_plus(diag,g,p,mu,eta,nm) result(tau)
      real(dp),intent(in)::diag(:),g(:),p(:),mu(:),eta;integer,intent(in)::nm
      real(dp)::a,b,c,rad,ts
      ts=diag(nm+1);a=ts+dot_r(diag,g,g,nm);b=dot_r(diag,mu,g,nm)-2.0_dp*dot_r(diag,p,g,nm)-eta*ts
      c=dot_r(diag,p,p,nm)-dot_r(diag,p,mu,nm);rad=b*b-4.0_dp*a*c;tau=(-b+sqrt(max(rad,0.0_dp)))/(2.0_dp*a)
   end function root_plus

   subroutine project_lin_sys(fac,g,diag,v,ut,n,m,iter)
      type(scs_ldlt_factor),intent(in)::fac;real(dp),intent(in)::g(:),diag(:),v(:);real(dp),intent(out)::ut(:)
      integer,intent(in)::n,m,iter;integer::nm
      nm=n+m;ut=v;ut(1:n)=ut(1:n)*diag(1:n);ut(n+1:nm)=-ut(n+1:nm)*diag(n+1:nm);call fac%solve(ut(1:nm))
      if(iter<feasible_iters)then;ut(nm+1)=1.0_dp;else;ut(nm+1)=root_plus(diag,g,ut,v,v(nm+1),nm);end if
      ut(1:nm)=ut(1:nm)-g*ut(nm+1)
   end subroutine project_lin_sys

   subroutine normalize_iterate(v)
      real(dp),intent(inout)::v(:);real(dp)::vn
      vn=norm_2(v);if(vn>0.0_dp)v=v*(sqrt(real(size(v),dp))*iterate_norm/vn)
   end subroutine normalize_iterate

   subroutine populate_residual(d,scal,have_scal,u,rsk,r,iter,n,m)
      type(scs_data),intent(in)::d;type(scs_scaling),intent(in)::scal;logical,intent(in)::have_scal
      real(dp),intent(in)::u(:),rsk(:);type(residual_state),intent(inout)::r;integer,intent(in)::iter,n,m
      real(dp)::nan,nm_ax_s,nm_px,nm_aty
      nan=ieee_value(1.0_dp,ieee_quiet_nan);r%last_iter=iter;r%tau=abs(u(n+m+1));r%kap=abs(rsk(n+m+1))
      r%x=u(1:n);r%y=u(n+1:n+m);r%s=rsk(n+1:n+m)
      if(have_scal)then
         r%x=r%x*(scal%E/scal%dual_scale);r%y=r%y*(scal%D/scal%primal_scale);r%s=r%s/(scal%D*scal%dual_scale)
      end if
      r%ax=0.0_dp;call accum_by_a(d%A,r%x,r%ax);r%ax_s=r%ax+r%s;r%ax_s_btau=r%ax_s-d%b*r%tau
      r%px=0.0_dp;if(d%has_p)call accum_by_p(d%P,r%x,r%px)
      r%aty=0.0_dp;call accum_by_atrans(d%A,r%y,r%aty);r%px_aty_ctau=r%px+r%aty+d%c*r%tau
      r%xt_p_x_tau=dot_product(r%px,r%x);r%bty_tau=dot_product(r%y,d%b);r%ctx_tau=dot_product(r%x,d%c)
      r%bty=safe_div_pos(r%bty_tau,r%tau);r%ctx=safe_div_pos(r%ctx_tau,r%tau);r%xt_p_x=safe_div_pos(r%xt_p_x_tau,r%tau*r%tau)
      r%gap=abs(r%xt_p_x+r%ctx+r%bty);r%pobj=r%xt_p_x/2.0_dp+r%ctx;r%dobj=-r%xt_p_x/2.0_dp-r%bty
      r%res_pri=safe_div_pos(norm_inf(r%ax_s_btau),r%tau);r%res_dual=safe_div_pos(norm_inf(r%px_aty_ctau),r%tau)
      r%res_unbdd_a=nan;r%res_unbdd_p=nan;r%res_infeas=nan
      if(r%ctx_tau < -infeas_negativity_tol)then
         nm_ax_s = norm_inf(r%ax_s)
         nm_px = norm_inf(r%px)
         r%res_unbdd_a = safe_div_pos(nm_ax_s,-r%ctx_tau)
         r%res_unbdd_p = safe_div_pos(nm_px,-r%ctx_tau)
      end if
      if(r%bty_tau < -infeas_negativity_tol)then;nm_aty=norm_inf(r%aty);r%res_infeas=safe_div_pos(nm_aty,-r%bty_tau);end if
   end subroutine populate_residual

   integer(i4) function convergence_status(r,nmb,nmc,s) result(status)
      type(residual_state),intent(in)::r;real(dp),intent(in)::nmb,nmc;type(scs_settings),intent(in)::s
      real(dp)::grl,prl,drl
      status=scs_unfinished
      if(r%tau>0.0_dp)then
         grl=max(abs(r%xt_p_x),abs(r%ctx),abs(r%bty))
         prl=max(nmb*r%tau,norm_inf(r%s),norm_inf(r%ax))/r%tau
         drl=max(nmc*r%tau,norm_inf(r%px),norm_inf(r%aty))/r%tau
         if(r%res_pri < s%eps_abs+s%eps_rel*prl .and. r%res_dual < s%eps_abs+s%eps_rel*drl .and. &
            r%gap < s%eps_abs+s%eps_rel*grl)then;status=scs_solved;return;end if
      end if
      if(.not.ieee_is_nan(r%res_unbdd_a) .and. .not.ieee_is_nan(r%res_unbdd_p))then
         if(r%res_unbdd_a<s%eps_infeas .and. r%res_unbdd_p<s%eps_infeas)then;status=scs_unbounded;return;end if
      end if
      if(.not.ieee_is_nan(r%res_infeas))then;if(r%res_infeas<s%eps_infeas)status=scs_infeasible;end if
   end function convergence_status

   subroutine maybe_update_scale(d,k,s,fac,diag,g,h,u,ut,v,rsk,r,nmb,nmc, &
      iter,last_update,sumlog,nlog,nupdates,aa)
      type(scs_data),intent(in)::d
      type(scs_cone),intent(in)::k
      type(scs_settings),intent(inout)::s;type(scs_ldlt_factor),intent(inout)::fac;real(dp),intent(inout)::diag(:),g(:),v(:)
      real(dp),intent(in)::h(:),u(:),ut(:),rsk(:),nmb,nmc;type(residual_state),intent(in)::r
      integer,intent(in)::iter;integer,intent(inout)::last_update,nlog,nupdates;real(dp),intent(inout)::sumlog
      type(aa_workspace),intent(inout)::aa
      real(dp)::denp,dend,relp,reld,factor,newscale
      logical::ok
      denp=max(norm_inf(r%ax),norm_inf(r%s),nmb*r%tau);relp=safe_div_pos(norm_inf(r%ax_s_btau),denp)
      dend=max(norm_inf(r%px),norm_inf(r%aty),nmc*r%tau);reld=safe_div_pos(norm_inf(r%px_aty_ctau),dend)
      if(relp<=0.0_dp .or. reld<=0.0_dp)return
      sumlog=sumlog+log(relp)-log(reld);nlog=nlog+1;factor=sqrt(exp(sumlog/real(nlog,dp)))
      if(iter-last_update<rescaling_min_iters)return
      newscale=min(max(s%scale*factor,min_scale_value),max_scale_value)
      if(abs(newscale-s%scale)<=epsilon(1.0_dp)*max(1.0_dp,abs(s%scale)))return
      if(factor>sqrt(10.0_dp) .or. factor<1.0_dp/sqrt(10.0_dp))then
         nupdates=nupdates+1;sumlog=0.0_dp;nlog=0;last_update=iter;s%scale=newscale
         call set_diag_r(k,s,diag,d%n,d%m);call fac%factorize(d%A,d%P,d%has_p,diag,ok);if(.not.ok)return
         call update_g(fac,h,g,d%n,d%m);v=rsk/diag+2.0_dp*ut-u
         if(s%acceleration_lookback/=0)call aa%reset()
      end if
   end subroutine maybe_update_scale

   subroutine finalize_solution(r,status_in,time_limit,sol,info)
      type(residual_state),intent(in)::r
      integer(i4),intent(in)::status_in
      logical,intent(in)::time_limit;type(scs_solution),intent(inout)::sol;type(scs_info),intent(inout)::info
      integer(i4)::status;real(dp)::nan,pinf,ninf
      nan=ieee_value(1.0_dp,ieee_quiet_nan);pinf=ieee_value(1.0_dp,ieee_positive_inf);ninf=ieee_value(1.0_dp,ieee_negative_inf)
      if(allocated(sol%x))deallocate(sol%x);if(allocated(sol%y))deallocate(sol%y);if(allocated(sol%s))deallocate(sol%s)
      allocate(sol%x(size(r%x)),sol%y(size(r%y)),sol%s(size(r%s)));sol%x=r%x;sol%y=r%y;sol%s=r%s
      status=status_in
      if(status==scs_unfinished)then
         if (r%tau > r%kap) then
            status = scs_solved_inaccurate
         else if (r%bty_tau < r%ctx_tau) then
            status = scs_infeasible_inaccurate
         else
            status = scs_unbounded_inaccurate
         end if
      end if
      select case(status)
      case(scs_solved,scs_solved_inaccurate)
         sol%x=sol%x*safe_div_pos(1.0_dp,r%tau);sol%y=sol%y*safe_div_pos(1.0_dp,r%tau);sol%s=sol%s*safe_div_pos(1.0_dp,r%tau)
         info%gap=r%gap;info%res_pri=r%res_pri;info%res_dual=r%res_dual;info%pobj=r%pobj;info%dobj=r%dobj;info%status='solved'
      case(scs_infeasible,scs_infeasible_inaccurate)
         sol%y = sol%y*(-1.0_dp/r%bty_tau)
         sol%x=nan;sol%s=nan;info%gap=nan;info%res_pri=nan;info%res_dual=nan
         info%pobj=pinf;info%dobj=pinf;info%status='infeasible'
      case(scs_unbounded,scs_unbounded_inaccurate)
         sol%x = sol%x*(-1.0_dp/r%ctx_tau)
         sol%s = sol%s*(-1.0_dp/r%ctx_tau)
         sol%y=nan;info%gap=nan;info%res_pri=nan;info%res_dual=nan
         info%pobj=ninf;info%dobj=ninf;info%status='unbounded'
      end select
      if(status==scs_solved_inaccurate .or. status==scs_infeasible_inaccurate .or. status==scs_unbounded_inaccurate)then
         if (time_limit) then
            info%status = trim(info%status)//' (inaccurate - reached time_limit_secs)'
         else
            info%status = trim(info%status)//' (inaccurate - reached max_iters)'
         end if
      end if
      info%status_val=status;info%res_infeas=r%res_infeas;info%res_unbdd_a=r%res_unbdd_a;info%res_unbdd_p=r%res_unbdd_p
   end subroutine finalize_solution
end module scs_solver
