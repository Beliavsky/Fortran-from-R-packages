! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_adjpin
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use pinstimation_kinds, only : dp, i8
   use pinstimation_types, only : trade_counts, adjpin_parameters, adjpin_restrictions, adjpin_result
   use pinstimation_math, only : logistic, logit, softplus, inv_softplus, log_sum_exp, poisson_log_pmf, random_poisson
   use pinstimation_optimization, only : optimizer_result, minimize_nelder_mead, minimize_bfgs
   implicit none
   private
   public :: adjpin_loglik, adjpin_distribution, adjpin_values
   public :: fit_adjpin, fit_adjpin_ml, fit_adjpin_ecm, simulate_adjpin, initial_adjpin

contains

   pure logical function valid_adjpin_parameters(p) result(ok)
      type(adjpin_parameters), intent(in) :: p
      ok = p%alpha >= 0.0_dp .and. p%alpha <= 1.0_dp .and. p%delta >= 0.0_dp .and. p%delta <= 1.0_dp .and. &
         p%theta >= 0.0_dp .and. p%theta <= 1.0_dp .and. p%theta_p >= 0.0_dp .and. p%theta_p <= 1.0_dp .and. &
         p%eps_b > 0.0_dp .and. p%eps_s > 0.0_dp .and. p%mu_b > 0.0_dp .and. p%mu_s > 0.0_dp .and. &
         p%d_b > 0.0_dp .and. p%d_s > 0.0_dp
   end function valid_adjpin_parameters

   pure function adjpin_distribution(p) result(w)
      type(adjpin_parameters), intent(in) :: p
      real(dp) :: w(6)
      w = [(1.0_dp-p%alpha)*(1.0_dp-p%theta), (1.0_dp-p%alpha)*p%theta, &
         p%alpha*(1.0_dp-p%delta)*(1.0_dp-p%theta_p), p%alpha*(1.0_dp-p%delta)*p%theta_p, &
         p%alpha*p%delta*(1.0_dp-p%theta_p), p%alpha*p%delta*p%theta_p]
   end function adjpin_distribution

   pure subroutine adjpin_values(p, adjpin, psos)
      type(adjpin_parameters), intent(in) :: p
      real(dp), intent(out) :: adjpin, psos
      real(dp) :: informed, shock, denominator
      informed = p%alpha*((1.0_dp-p%delta)*p%mu_b + p%delta*p%mu_s)
      shock = (p%d_b+p%d_s)*(p%alpha*p%theta_p + (1.0_dp-p%alpha)*p%theta)
      denominator = informed + shock + p%eps_b + p%eps_s
      adjpin = informed/denominator
      psos = shock/denominator
   end subroutine adjpin_values

   real(dp) function adjpin_loglik(data, p) result(value)
      type(trade_counts), intent(in) :: data
      type(adjpin_parameters), intent(in) :: p
      real(dp) :: w(6), terms(6)
      integer :: i
      if (.not. data%valid() .or. .not. valid_adjpin_parameters(p)) then
         value = -huge(1.0_dp)
         return
      end if
      w = adjpin_distribution(p)
      value = 0.0_dp
      do i = 1, data%size()
         terms(1) = safe_log(w(1)) + poisson_log_pmf(data%buys(i),p%eps_b) + poisson_log_pmf(data%sells(i),p%eps_s)
         terms(2) = safe_log(w(2)) + poisson_log_pmf(data%buys(i),p%eps_b+p%d_b) + &
            poisson_log_pmf(data%sells(i),p%eps_s+p%d_s)
         terms(3) = safe_log(w(3)) + poisson_log_pmf(data%buys(i),p%eps_b+p%mu_b) + poisson_log_pmf(data%sells(i),p%eps_s)
         terms(4) = safe_log(w(4)) + poisson_log_pmf(data%buys(i),p%eps_b+p%d_b+p%mu_b) + &
            poisson_log_pmf(data%sells(i),p%eps_s+p%d_s)
         terms(5) = safe_log(w(5)) + poisson_log_pmf(data%buys(i),p%eps_b) + poisson_log_pmf(data%sells(i),p%eps_s+p%mu_s)
         terms(6) = safe_log(w(6)) + poisson_log_pmf(data%buys(i),p%eps_b+p%d_b) + &
            poisson_log_pmf(data%sells(i),p%eps_s+p%d_s+p%mu_s)
         value = value + log_sum_exp(terms)
      end do
   end function adjpin_loglik

   function initial_adjpin(data, restrictions) result(p)
      type(trade_counts), intent(in) :: data
      type(adjpin_restrictions), intent(in), optional :: restrictions
      type(adjpin_parameters) :: p
      type(adjpin_restrictions) :: r
      real(dp), allocatable :: b(:), s(:)
      real(dp) :: mb, ms, imbalance
      integer :: n
      r = adjpin_restrictions()
      if (present(restrictions)) r = restrictions
      n = data%size()
      allocate(b(n),s(n))
      b = real(data%buys,dp); s = real(data%sells,dp)
      mb = sum(b)/real(max(1,n),dp); ms = sum(s)/real(max(1,n),dp)
      imbalance = sum(abs(b-s))/real(max(1,n),dp)
      p%alpha = 0.35_dp
      p%delta = min(0.8_dp,max(0.2_dp,real(count(s>b),dp)/real(max(1,n),dp)))
      p%theta = 0.15_dp; p%theta_p = 0.25_dp
      p%eps_b = max(0.1_dp,0.55_dp*mb); p%eps_s = max(0.1_dp,0.55_dp*ms)
      p%mu_b = max(0.2_dp,imbalance); p%mu_s = max(0.2_dp,imbalance)
      p%d_b = max(0.2_dp,0.2_dp*mb); p%d_s = max(0.2_dp,0.2_dp*ms)
      call apply_restrictions(p,r)
   end function initial_adjpin

   subroutine fit_adjpin(data, result, restrictions, initial, method, optimizer, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      type(adjpin_result), intent(out) :: result
      type(adjpin_restrictions), intent(in), optional :: restrictions
      type(adjpin_parameters), intent(in), optional :: initial
      character(len=*), intent(in), optional :: method, optimizer
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      character(len=16) :: selected
      selected = 'ML'
      if (present(method)) selected = uppercase(trim(method))
      if (trim(selected) == 'ECM') then
         call fit_adjpin_ecm(data,result,restrictions,initial,max_iterations,tolerance)
      else
         call fit_adjpin_ml(data,result,restrictions,initial,optimizer,max_iterations,tolerance)
      end if
   end subroutine fit_adjpin

   subroutine fit_adjpin_ml(data, result, restrictions, initial, method, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      type(adjpin_result), intent(out) :: result
      type(adjpin_restrictions), intent(in), optional :: restrictions
      type(adjpin_parameters), intent(in), optional :: initial
      character(len=*), intent(in), optional :: method
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(adjpin_restrictions) :: r
      type(adjpin_parameters) :: p0, p, starts(5)
      type(optimizer_result) :: opt, best
      real(dp), allocatable :: u0(:)
      character(len=16) :: solver
      integer :: i,maxit
      real(dp) :: tol
      r = adjpin_restrictions(); if (present(restrictions)) r = restrictions
      p0 = initial_adjpin(data,r); if (present(initial)) p0 = initial
      call apply_restrictions(p0,r)
      starts(1)=p0
      starts(2)=p0; starts(2)%alpha=0.15_dp; starts(2)%theta=0.35_dp; starts(2)%theta_p=0.35_dp
      starts(3)=p0; starts(3)%alpha=0.65_dp; starts(3)%delta=0.25_dp
      starts(4)=p0; starts(4)%alpha=0.65_dp; starts(4)%delta=0.75_dp
      starts(5)=p0; starts(5)%d_b=2.0_dp*p0%d_b; starts(5)%d_s=2.0_dp*p0%d_s
      do i=1,size(starts); call apply_restrictions(starts(i),r); end do
      solver='NELDER-MEAD'; if(present(method)) solver=uppercase(trim(method))
      maxit=5000; tol=1.0e-8_dp
      if(present(max_iterations)) maxit=max_iterations
      if(present(tolerance)) tol=tolerance
      best%objective=huge(1.0_dp)
      do i=1,size(starts)
         call adjpin_pack(starts(i),r,u0)
         if(solver(1:min(4,len_trim(solver)))=='BFGS') then
            call minimize_bfgs(objective,u0,opt,tol,maxit)
         else
            call minimize_nelder_mead(objective,u0,opt,tol,maxit)
         end if
         if(opt%objective<best%objective) best=opt
      end do
      call adjpin_unpack(best%parameters,r,p)
      result%parameters=p; result%log_likelihood=-best%objective
      result%iterations=best%iterations; result%evaluations=best%evaluations
      result%status=best%status; result%converged=best%converged
      call adjpin_values(p,result%adjpin,result%psos)
   contains
      real(dp) function objective(u) result(v)
         real(dp),intent(in)::u(:)
         type(adjpin_parameters)::pp
         call adjpin_unpack(u,r,pp)
         v=-adjpin_loglik(data,pp)
         if(.not.ieee_is_finite(v)) v=huge(1.0_dp)/100.0_dp
      end function objective
   end subroutine fit_adjpin_ml

   subroutine fit_adjpin_ecm(data, result, restrictions, initial, max_iterations, tolerance)
      type(trade_counts), intent(in) :: data
      type(adjpin_result), intent(out) :: result
      type(adjpin_restrictions), intent(in), optional :: restrictions
      type(adjpin_parameters), intent(in), optional :: initial
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(adjpin_restrictions) :: r
      type(adjpin_parameters) :: p
      type(optimizer_result) :: opt
      real(dp),allocatable::post(:,:),u(:)
      real(dp)::terms(6),den,oldll,newll,tol,w(6)
      integer::i,iter,maxit
      r=adjpin_restrictions(); if(present(restrictions)) r=restrictions
      p=initial_adjpin(data,r); if(present(initial)) p=initial
      call apply_restrictions(p,r)
      maxit=120; tol=1.0e-6_dp
      if(present(max_iterations)) maxit=max_iterations
      if(present(tolerance)) tol=tolerance
      allocate(post(data%size(),6)); oldll=adjpin_loglik(data,p)
      do iter=1,maxit
         w=adjpin_distribution(p)
         do i=1,data%size()
            terms(1)=safe_log(w(1))+poisson_log_pmf(data%buys(i),p%eps_b)+poisson_log_pmf(data%sells(i),p%eps_s)
            terms(2)=safe_log(w(2))+poisson_log_pmf(data%buys(i),p%eps_b+p%d_b)+poisson_log_pmf(data%sells(i),p%eps_s+p%d_s)
            terms(3)=safe_log(w(3))+poisson_log_pmf(data%buys(i),p%eps_b+p%mu_b)+poisson_log_pmf(data%sells(i),p%eps_s)
            terms(4)=safe_log(w(4))+poisson_log_pmf(data%buys(i),p%eps_b+p%d_b+p%mu_b)+poisson_log_pmf(data%sells(i),p%eps_s+p%d_s)
            terms(5)=safe_log(w(5))+poisson_log_pmf(data%buys(i),p%eps_b)+poisson_log_pmf(data%sells(i),p%eps_s+p%mu_s)
            terms(6)=safe_log(w(6))+poisson_log_pmf(data%buys(i),p%eps_b+p%d_b)+poisson_log_pmf(data%sells(i),p%eps_s+p%d_s+p%mu_s)
            den=log_sum_exp(terms); post(i,:)=exp(terms-den)
         end do
         w=sum(post,dim=1)/real(data%size(),dp)
         p%alpha=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,sum(w(3:6))))
         p%delta=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,(w(5)+w(6))/max(sum(w(3:6)),1.0e-12_dp)))
         p%theta=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,w(2)/max(w(1)+w(2),1.0e-12_dp)))
         p%theta_p=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,(w(4)+w(6))/max(sum(w(3:6)),1.0e-12_dp)))
         call apply_restrictions(p,r)
         call adjpin_pack(p,r,u)
         call minimize_nelder_mead(qobjective,u,opt,1.0e-6_dp,900)
         call adjpin_unpack(opt%parameters,r,p)
         newll=adjpin_loglik(data,p)
         if(abs(newll-oldll)<=tol*(1.0_dp+abs(oldll))) exit
         oldll=newll
      end do
      result%parameters=p; result%log_likelihood=adjpin_loglik(data,p)
      result%iterations=iter; result%evaluations=opt%evaluations
      result%converged=iter<=maxit; result%status=merge(0,1,result%converged)
      call adjpin_values(p,result%adjpin,result%psos)
   contains
      real(dp) function qobjective(uvec) result(v)
         real(dp),intent(in)::uvec(:)
         type(adjpin_parameters)::pp
         real(dp)::q
         integer::ii
         call adjpin_unpack(uvec,r,pp); q=0.0_dp
         w=adjpin_distribution(pp)
         do ii=1,data%size()
            q=q+post(ii,1)*(safe_log(w(1))+poisson_log_pmf(data%buys(ii),pp%eps_b)+poisson_log_pmf(data%sells(ii),pp%eps_s))
            q=q+post(ii,2)*(safe_log(w(2))+poisson_log_pmf(data%buys(ii),pp%eps_b+pp%d_b)+poisson_log_pmf(data%sells(ii),pp%eps_s+pp%d_s))
            q=q+post(ii,3)*(safe_log(w(3))+poisson_log_pmf(data%buys(ii),pp%eps_b+pp%mu_b)+poisson_log_pmf(data%sells(ii),pp%eps_s))
            q=q+post(ii,4)*(safe_log(w(4))+poisson_log_pmf(data%buys(ii),pp%eps_b+pp%d_b+pp%mu_b)+poisson_log_pmf(data%sells(ii),pp%eps_s+pp%d_s))
            q=q+post(ii,5)*(safe_log(w(5))+poisson_log_pmf(data%buys(ii),pp%eps_b)+poisson_log_pmf(data%sells(ii),pp%eps_s+pp%mu_s))
            q=q+post(ii,6)*(safe_log(w(6))+poisson_log_pmf(data%buys(ii),pp%eps_b+pp%d_b)+poisson_log_pmf(data%sells(ii),pp%eps_s+pp%d_s+pp%mu_s))
         end do
         v=-q
         if(.not.ieee_is_finite(v)) v=huge(1.0_dp)/100.0_dp
      end function qobjective
   end subroutine fit_adjpin_ecm

   subroutine simulate_adjpin(days,p,data,states,seed,status)
      integer,intent(in)::days
      type(adjpin_parameters),intent(in)::p
      type(trade_counts),intent(out)::data
      integer,allocatable,intent(out),optional::states(:)
      integer,intent(in),optional::seed
      integer,intent(out),optional::status
      real(dp)::w(6),cw(6),u,lb,ls
      integer::i,state,st1,st2
      if(present(status)) status=0
      if(days<1.or..not.valid_adjpin_parameters(p)) then
         allocate(data%buys(0),data%sells(0)); if(present(status)) status=1; return
      end if
      if(present(seed)) call seed_local(seed)
      allocate(data%buys(days),data%sells(days)); if(present(states)) allocate(states(days))
      w=adjpin_distribution(p); cw(1)=w(1)
      do i=2,6; cw(i)=cw(i-1)+w(i); end do
      do i=1,days
         call random_number(u); state=1
         do while(state<6.and.u>cw(state)); state=state+1; end do
         select case(state)
         case(1); lb=p%eps_b; ls=p%eps_s
         case(2); lb=p%eps_b+p%d_b; ls=p%eps_s+p%d_s
         case(3); lb=p%eps_b+p%mu_b; ls=p%eps_s
         case(4); lb=p%eps_b+p%d_b+p%mu_b; ls=p%eps_s+p%d_s
         case(5); lb=p%eps_b; ls=p%eps_s+p%mu_s
         case default; lb=p%eps_b+p%d_b; ls=p%eps_s+p%d_s+p%mu_s
         end select
         data%buys(i)=random_poisson(lb,st1); data%sells(i)=random_poisson(ls,st2)
         if(present(states)) states(i)=state
         if(present(status)) status=max(status,max(st1,st2))
      end do
   end subroutine simulate_adjpin

   subroutine adjpin_pack(p,r,u)
      type(adjpin_parameters),intent(in)::p
      type(adjpin_restrictions),intent(in)::r
      real(dp),allocatable,intent(out)::u(:)
      integer::n,k
      n=10-merge(1,0,r%equal_theta)-merge(1,0,r%equal_eps)-merge(1,0,r%equal_mu)-merge(1,0,r%equal_d)
      allocate(u(n)); k=0
      k=k+1; u(k)=logit(p%alpha); k=k+1; u(k)=logit(p%delta); k=k+1; u(k)=logit(p%theta)
      if(.not.r%equal_theta) then; k=k+1; u(k)=logit(p%theta_p); end if
      k=k+1; u(k)=inv_softplus(p%eps_b); if(.not.r%equal_eps) then; k=k+1; u(k)=inv_softplus(p%eps_s); end if
      k=k+1; u(k)=inv_softplus(p%mu_b); if(.not.r%equal_mu) then; k=k+1; u(k)=inv_softplus(p%mu_s); end if
      k=k+1; u(k)=inv_softplus(p%d_b); if(.not.r%equal_d) then; k=k+1; u(k)=inv_softplus(p%d_s); end if
   end subroutine adjpin_pack

   subroutine adjpin_unpack(u,r,p)
      real(dp),intent(in)::u(:)
      type(adjpin_restrictions),intent(in)::r
      type(adjpin_parameters),intent(out)::p
      integer::k
      k=0; k=k+1;p%alpha=logistic(u(k)); k=k+1;p%delta=logistic(u(k)); k=k+1;p%theta=logistic(u(k))
      if(r%equal_theta) then;p%theta_p=p%theta;else;k=k+1;p%theta_p=logistic(u(k));end if
      k=k+1;p%eps_b=softplus(u(k))+1.0e-10_dp
      if(r%equal_eps) then;p%eps_s=p%eps_b;else;k=k+1;p%eps_s=softplus(u(k))+1.0e-10_dp;end if
      k=k+1;p%mu_b=softplus(u(k))+1.0e-10_dp
      if(r%equal_mu) then;p%mu_s=p%mu_b;else;k=k+1;p%mu_s=softplus(u(k))+1.0e-10_dp;end if
      k=k+1;p%d_b=softplus(u(k))+1.0e-10_dp
      if(r%equal_d) then;p%d_s=p%d_b;else;k=k+1;p%d_s=softplus(u(k))+1.0e-10_dp;end if
   end subroutine adjpin_unpack

   pure subroutine apply_restrictions(p,r)
      type(adjpin_parameters),intent(inout)::p
      type(adjpin_restrictions),intent(in)::r
      if(r%equal_theta) p%theta_p=p%theta
      if(r%equal_eps) p%eps_s=p%eps_b
      if(r%equal_mu) p%mu_s=p%mu_b
      if(r%equal_d) p%d_s=p%d_b
   end subroutine apply_restrictions

   pure real(dp) function safe_log(x) result(v)
      real(dp),intent(in)::x
      if(x>0.0_dp) then;v=log(x);else;v=-huge(1.0_dp);end if
   end function safe_log

   pure function uppercase(text) result(out)
      character(len=*),intent(in)::text
      character(len=len(text))::out
      integer::i,code
      do i=1,len(text);code=iachar(text(i:i));if(code>=iachar('a').and.code<=iachar('z')) then;out(i:i)=achar(code-32);else;out(i:i)=text(i:i);end if;end do
   end function uppercase

   subroutine seed_local(seed)
      integer,intent(in)::seed
      integer,allocatable::put(:)
      integer::n,i
      call random_seed(size=n);allocate(put(n))
      do i=1,n;put(i)=modulo(seed+15485863*i,huge(1)-1);if(put(i)<=0)put(i)=i;end do
      call random_seed(put=put)
   end subroutine seed_local

end module pinstimation_adjpin
