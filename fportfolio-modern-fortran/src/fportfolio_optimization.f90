! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_optimization
  use fportfolio_kinds, only: dp
  use fportfolio_types, only: linear_constraints, optimizer_result, frontier_result
  use fportfolio_linalg, only: quadratic_form, symmetric_eigen
  use fportfolio_risk, only: covariance_risk, historical_es, portfolio_returns, diversification_ratio
  implicit none
  private
  public :: initialize_constraints, feasible_portfolio, solve_quadratic_program, solve_linear_program, &
            minvariance_portfolio, efficient_portfolio, maxreturn_portfolio, tangency_portfolio, &
            maxratio_portfolio, minimum_mad_portfolio, minimum_cvar_portfolio, risk_parity_portfolio, &
            maximum_diversification_portfolio, portfolio_frontier, cardinality_minvariance_portfolio, &
            project_feasible
contains
  subroutine initialize_constraints(n,con,lower,upper,budget)
    integer,intent(in)::n
    type(linear_constraints),intent(out)::con
    real(dp),intent(in),optional::lower(:),upper(:),budget
    allocate(con%lower(n),con%upper(n))
    con%lower=0.0_dp;con%upper=1.0_dp
    if(present(lower))con%lower=lower
    if(present(upper))con%upper=upper
    if(present(budget))then;con%budget=budget;con%budget_min=budget;con%budget_max=budget;end if
  end subroutine initialize_constraints

  subroutine feasible_portfolio(con,weights,status)
    type(linear_constraints),intent(in)::con
    real(dp),allocatable,intent(out)::weights(:)
    integer,intent(out)::status
    real(dp),allocatable::y(:)
    integer::n
    n=size(con%lower);allocate(y(n),weights(n));y=0.5_dp*(con%lower+con%upper)
    call project_feasible(y,con,weights,status)
  end subroutine feasible_portfolio

  subroutine project_feasible(y,con,x,status,mu)
    real(dp),intent(in)::y(:)
    type(linear_constraints),intent(in)::con
    real(dp),allocatable,intent(out)::x(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::mu(:)
    real(dp),allocatable::z(:),old(:),a(:)
    integer::iter,n,j
    real(dp)::den,diff,target
    n=size(y)
    allocate(x(n),z(n),old(n))
    x=y
    status=0
    do iter=1,1000
      old=x
      if(con%budget_equality) then
        call project_box_sum(x,con%lower,con%upper,con%budget,z,status)
        if(status/=0)return
        x=z
      else
        x=min(con%upper,max(con%lower,x))
        if(sum(x)<con%budget_min)then
          call project_box_sum(x,con%lower,con%upper,con%budget_min,z,status);x=z
        else if(sum(x)>con%budget_max)then
          call project_box_sum(x,con%lower,con%upper,con%budget_max,z,status);x=z
        end if
      end if
      if(con%has_target_return .and. present(mu))then
        a=mu;den=dot_product(a,a)
        if(den>0.0_dp)x=x+a*(con%target_return-dot_product(a,x))/den
        x=min(con%upper,max(con%lower,x))
      end if
      if(allocated(con%a_eq))then
        do j=1,size(con%a_eq,1)
          den=dot_product(con%a_eq(j,:),con%a_eq(j,:))
          if(den>0.0_dp)x=x+con%a_eq(j,:)*(con%b_eq(j)-dot_product(con%a_eq(j,:),x))/den
          x=min(con%upper,max(con%lower,x))
        end do
      end if
      if(allocated(con%a_ineq))then
        do j=1,size(con%a_ineq,1)
          diff=dot_product(con%a_ineq(j,:),x)-con%b_ineq(j)
          if(diff>0.0_dp)then
            den=dot_product(con%a_ineq(j,:),con%a_ineq(j,:))
            if(den>0.0_dp)x=x-con%a_ineq(j,:)*diff/den
            x=min(con%upper,max(con%lower,x))
          end if
        end do
      end if
      if(maxval(abs(x-old))<1.0e-12_dp)exit
    end do
    if(con%budget_equality .and. abs(sum(x)-con%budget)>1.0e-7_dp)status=2
    if(con%has_target_return .and. present(mu))then
      target=dot_product(mu,x)
      if(abs(target-con%target_return)>1.0e-6_dp)status=3
    end if
  end subroutine project_feasible

  subroutine project_box_sum(y,lower,upper,total,x,status)
    real(dp),intent(in)::y(:),lower(:),upper(:),total
    real(dp),allocatable,intent(out)::x(:)
    integer,intent(out)::status
    real(dp)::lo,hi,mid,s
    integer::it,n
    n=size(y);allocate(x(n));status=0
    if(total<sum(lower)-1.0e-12_dp .or. total>sum(upper)+1.0e-12_dp)then
      x=min(upper,max(lower,y));status=1;return
    end if
    lo=minval(y-upper)-1.0_dp;hi=maxval(y-lower)+1.0_dp
    do it=1,200
      mid=0.5_dp*(lo+hi)
      x=min(upper,max(lower,y-mid))
      s=sum(x)
      if(abs(s-total)<1.0e-13_dp)return
      if(s>total)then
        lo=mid
      else
        hi=mid
      end if
    end do
    x=min(upper,max(lower,y-0.5_dp*(lo+hi)))
  end subroutine project_box_sum

  subroutine solve_quadratic_program(q,c,con,result,mu,max_iter,tol)
    real(dp),intent(in)::q(:,:),c(:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    real(dp),intent(in),optional::mu(:),tol
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),grad(:),vals(:),vecs(:,:),y0(:)
    real(dp)::step,tolerance,obj,objn,lipschitz
    integer::iter,n,status,mi,info
    n=size(c);mi=5000;if(present(max_iter))mi=max_iter;tolerance=1.0e-10_dp;if(present(tol))tolerance=tol
    if(present(mu) .and. con%has_target_return)then
      allocate(y0(n))
      y0=0.5_dp*(con%lower+con%upper)
      call project_feasible(y0,con,w,status,mu)
    else
      call feasible_portfolio(con,w,status)
    end if
    if(status/=0)then
      call set_failure(result,n,"infeasible constraints",status)
      return
    end if
    call symmetric_eigen(q,vals,vecs,info)
    if(info==0)then;lipschitz=max(maxval(abs(vals)),1.0e-8_dp);else;lipschitz=max(maxval(abs(q)),1.0_dp)*real(n,dp);end if
    step=1.0_dp/lipschitz;allocate(grad(n));obj=0.5_dp*quadratic_form(w,q)+dot_product(c,w)
    do iter=1,mi
      grad=matmul(q,w)+c
      if(present(mu))then;call project_feasible(w-step*grad,con,wn,status,mu);else;call project_feasible(w-step*grad,con,wn,status);end if
      if(status/=0)then;step=0.5_dp*step;cycle;end if
      objn=0.5_dp*quadratic_form(wn,q)+dot_product(c,wn)
      if(objn>obj+1.0e-12_dp)then;step=0.5_dp*step;cycle;end if
      if(maxval(abs(wn-w))<tolerance)exit
      w=wn;obj=objn;step=min(step*1.02_dp,1.0_dp/lipschitz)
    end do
    call fill_result(result,w,obj,iter,iter,iter<mi,"quadratic program")
    if(present(mu))result%expected_return=dot_product(mu,w)
    result%risk=sqrt(max(quadratic_form(w,q),0.0_dp))
  end subroutine solve_quadratic_program

  subroutine solve_linear_program(c,con,result,mu,max_iter,tol)
    real(dp),intent(in)::c(:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    real(dp),intent(in),optional::mu(:),tol
    integer,intent(in),optional::max_iter
    real(dp),allocatable::q(:,:)
    integer::n
    n=size(c);allocate(q(n,n));q=0.0_dp
    call solve_quadratic_program(q,c,con,result,mu,max_iter,tol)
    result%message="linear program by projected active constraints"
  end subroutine solve_linear_program

  subroutine minvariance_portfolio(mu,sigma,con,result)
    real(dp),intent(in)::mu(:),sigma(:,:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    real(dp)::c(size(mu));c=0.0_dp
    call solve_quadratic_program(sigma,c,con,result,mu=mu)
    result%objective=result%risk
    result%message="minimum variance portfolio"
  end subroutine minvariance_portfolio

  subroutine efficient_portfolio(mu,sigma,target_return,con,result)
    real(dp),intent(in)::mu(:),sigma(:,:),target_return
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    type(linear_constraints)::ct
    real(dp)::c(size(mu));c=0.0_dp;ct=con;ct%has_target_return=.true.;ct%target_return=target_return
    call solve_quadratic_program(sigma,c,ct,result,mu=mu)
    result%objective=result%risk;result%message="target-return efficient portfolio"
  end subroutine efficient_portfolio

  subroutine maxreturn_portfolio(mu,sigma,con,result)
    real(dp),intent(in)::mu(:),sigma(:,:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    call solve_linear_program(-mu,con,result,mu=mu)
    result%objective=-result%expected_return;result%risk=covariance_risk(sigma,result%weights)
    result%message="maximum return portfolio"
  end subroutine maxreturn_portfolio

  subroutine tangency_portfolio(mu,sigma,risk_free,con,result,max_iter)
    real(dp),intent(in)::mu(:),sigma(:,:),risk_free
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),grad(:),excess(:)
    real(dp)::s,a,step,sh,shn
    integer::it,mi,status,n
    n=size(mu);mi=10000;if(present(max_iter))mi=max_iter;allocate(excess(n),grad(n));excess=mu-risk_free
    call feasible_portfolio(con,w,status)
    if(status/=0)then
      call set_failure(result,n,"infeasible constraints",status)
      return
    end if
    step=0.05_dp
    sh=-huge(1.0_dp)
    do it=1,mi
      s=covariance_risk(sigma,w);a=dot_product(excess,w)
      if(s<=1.0e-14_dp)exit
      grad=excess/s-a*matmul(sigma,w)/(s**3)
      call project_feasible(w+step*grad,con,wn,status)
      if(status/=0)then;step=0.5_dp*step;cycle;end if
      shn=dot_product(excess,wn)/max(covariance_risk(sigma,wn),1.0e-14_dp)
      if(shn<sh)then;step=0.5_dp*step;cycle;end if
      if(maxval(abs(wn-w))<1.0e-11_dp)exit
      w=wn;sh=shn;step=min(0.2_dp,step*1.005_dp)
    end do
    call fill_result(result,w,-sh,it,it,it<mi,"tangency portfolio")
    result%expected_return=dot_product(mu,w);result%risk=covariance_risk(sigma,w)
    result%sharpe=(result%expected_return-risk_free)/max(result%risk,1.0e-14_dp)
  end subroutine tangency_portfolio

  subroutine maxratio_portfolio(mu,sigma,risk_free,con,result)
    real(dp),intent(in)::mu(:),sigma(:,:),risk_free
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    call tangency_portfolio(mu,sigma,risk_free,con,result)
  end subroutine maxratio_portfolio

  subroutine maximum_diversification_portfolio(mu,sigma,con,result,max_iter)
    real(dp),intent(in)::mu(:),sigma(:,:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),grad(:),sd(:)
    real(dp)::pvol,num,dr,drn,step
    integer::n,it,mi,status,i
    n=size(mu);mi=10000;if(present(max_iter))mi=max_iter;allocate(grad(n),sd(n))
    do i=1,n;sd(i)=sqrt(max(sigma(i,i),0.0_dp));end do
    call feasible_portfolio(con,w,status);step=0.05_dp;dr=diversification_ratio(sigma,w)
    do it=1,mi
      pvol=covariance_risk(sigma,w);num=dot_product(sd,w)
      if(pvol<=1.0e-14_dp)exit
      grad=sd/pvol-num*matmul(sigma,w)/(pvol**3)
      call project_feasible(w+step*grad,con,wn,status)
      drn=diversification_ratio(sigma,wn)
      if(drn<dr)then;step=0.5_dp*step;cycle;end if
      if(maxval(abs(wn-w))<1.0e-11_dp)exit
      w=wn;dr=drn;step=min(0.2_dp,step*1.005_dp)
    end do
    call fill_result(result,w,-dr,it,it,it<mi,"maximum diversification portfolio")
    result%expected_return=dot_product(mu,w);result%risk=covariance_risk(sigma,w)
  end subroutine maximum_diversification_portfolio

  subroutine risk_parity_portfolio(mu,sigma,con,result,target_budget,max_iter)
    real(dp),intent(in)::mu(:),sigma(:,:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    real(dp),intent(in),optional::target_budget(:)
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),b(:),m(:),rc(:)
    real(dp)::risk,err,eta
    integer::n,it,mi,status
    n=size(mu);mi=20000;if(present(max_iter))mi=max_iter;allocate(b(n),m(n),rc(n))
    if(present(target_budget))then;b=target_budget/sum(target_budget);else;b=1.0_dp/real(n,dp);end if
    call feasible_portfolio(con,w,status);w=max(w,1.0e-10_dp);call project_feasible(w,con,wn,status);w=wn
    eta=0.25_dp
    do it=1,mi
      risk=covariance_risk(sigma,w);m=matmul(sigma,w)/max(risk,1.0e-14_dp);rc=w*m
      err=maxval(abs(rc-b*risk))
      if(err<1.0e-10_dp)exit
      w=w*exp(max(-50.0_dp,min(50.0_dp,-eta*(rc/max(b*risk,1.0e-14_dp)-1.0_dp))))
      call project_feasible(w,con,wn,status);w=max(wn,1.0e-14_dp)
      eta=max(0.01_dp,eta*0.9995_dp)
    end do
    risk=covariance_risk(sigma,w);m=matmul(sigma,w)/max(risk,1.0e-14_dp);rc=w*m
    call fill_result(result,w,sum((rc-b*risk)**2),it,it,it<mi,"risk parity portfolio")
    result%expected_return=dot_product(mu,w);result%risk=risk
  end subroutine risk_parity_portfolio

  subroutine minimum_mad_portfolio(data,mu,con,result,max_iter)
    real(dp),intent(in)::data(:,:),mu(:)
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),r(:),rn(:),grad(:),centered(:,:)
    real(dp)::obj,objn,step,eps_smooth
    integer::n,it,mi,status
    n=size(mu)
    mi=20000
    if(present(max_iter))mi=max_iter
    allocate(grad(n),r(size(data,1)),rn(size(data,1)),centered(size(data,1),n))
    centered=data-spread(mu,1,size(data,1))
    call feasible_portfolio(con,w,status)
    step=0.5_dp
    eps_smooth=1.0e-6_dp
    r=matmul(centered,w)
    obj=sum(sqrt(r*r+eps_smooth**2))/real(size(r),dp)
    do it=1,mi
      grad=matmul(transpose(centered),r/sqrt(r*r+eps_smooth**2))/real(size(r),dp)
      call project_feasible(w-step*grad,con,wn,status,mu)
      rn=matmul(centered,wn)
      objn=sum(sqrt(rn*rn+eps_smooth**2))/real(size(rn),dp)
      if(objn>obj)then
        step=0.5_dp*step
        if(step<1.0e-12_dp)exit
        cycle
      end if
      if(maxval(abs(wn-w))<1.0e-10_dp .or. abs(objn-obj)<1.0e-12_dp)then
        w=wn
        obj=objn
        exit
      end if
      w=wn
      r=rn
      obj=objn
      step=min(1.0_dp,step*1.02_dp)
    end do
    r=matmul(centered,w)
    obj=sum(abs(r))/real(size(r),dp)
    call fill_result(result,w,obj,it,it,it<=mi,"minimum MAD portfolio")
    result%expected_return=dot_product(mu,w)
    result%risk=obj
  end subroutine minimum_mad_portfolio

  subroutine minimum_cvar_portfolio(data,mu,alpha,con,result,max_iter)
    real(dp),intent(in)::data(:,:),mu(:),alpha
    type(linear_constraints),intent(in)::con
    type(optimizer_result),intent(out)::result
    integer,intent(in),optional::max_iter
    real(dp),allocatable::w(:),wn(:),wbar(:),q(:),qn(:),lower_q(:),upper_q(:)
    real(dp)::tau,sigma_step,norm_x,cap,obj,change
    integer::n,p,it,mi,status
    n=size(data,1)
    p=size(data,2)
    mi=100000
    if(present(max_iter))mi=max_iter
    allocate(wbar(p),q(n),lower_q(n),upper_q(n))
    call feasible_portfolio(con,w,status)
    if(status/=0)then
      call set_failure(result,p,"infeasible constraints",status)
      return
    end if
    wbar=w
    cap=1.0_dp/max(alpha*real(n,dp),1.0_dp)
    lower_q=0.0_dp
    upper_q=cap
    q=1.0_dp/real(n,dp)
    norm_x=sqrt(sum(data*data))
    if(norm_x<=tiny(1.0_dp))then
      call fill_result(result,w,0.0_dp,0,0,.true.,"minimum CVaR portfolio")
      result%expected_return=dot_product(mu,w)
      result%risk=0.0_dp
      return
    end if
    tau=0.99_dp/norm_x
    sigma_step=0.99_dp/norm_x
    do it=1,mi
      call project_box_sum(q-sigma_step*matmul(data,wbar),lower_q,upper_q,1.0_dp,qn,status)
      if(status/=0)exit
      call project_feasible(w+tau*matmul(transpose(data),qn),con,wn,status,mu)
      if(status/=0)exit
      change=max(maxval(abs(wn-w)),maxval(abs(qn-q)))
      wbar=2.0_dp*wn-w
      w=wn
      q=qn
      if(change<1.0e-8_dp)exit
    end do
    obj=historical_es(data,w,alpha)
    call fill_result(result,w,obj,it,it,it<=mi .and. status==0,"minimum CVaR portfolio")
    result%expected_return=dot_product(mu,w)
    result%risk=obj
  end subroutine minimum_cvar_portfolio

  subroutine portfolio_frontier(mu,sigma,con,npoints,frontier)
    real(dp),intent(in)::mu(:),sigma(:,:)
    type(linear_constraints),intent(in)::con
    integer,intent(in)::npoints
    type(frontier_result),intent(out)::frontier
    type(optimizer_result)::lo,hi,res
    real(dp)::rmin,rmax
    integer::i,n
    n=size(mu)
    allocate(frontier%target_return(npoints),frontier%risk(npoints))
    allocate(frontier%weights(n,npoints),frontier%feasible(npoints))
    call minvariance_portfolio(mu,sigma,con,lo);call maxreturn_portfolio(mu,sigma,con,hi)
    rmin=lo%expected_return;rmax=hi%expected_return
    do i=1,npoints
      frontier%target_return(i)=rmin+(rmax-rmin)*real(i-1,dp)/real(max(1,npoints-1),dp)
      call efficient_portfolio(mu,sigma,frontier%target_return(i),con,res)
      frontier%feasible(i)=res%converged
      if(allocated(res%weights))frontier%weights(:,i)=res%weights
      frontier%risk(i)=res%risk
    end do
  end subroutine portfolio_frontier

  subroutine cardinality_minvariance_portfolio(mu,sigma,con,max_assets,min_buyin,result)
    real(dp),intent(in)::mu(:),sigma(:,:),min_buyin
    type(linear_constraints),intent(in)::con
    integer,intent(in)::max_assets
    type(optimizer_result),intent(out)::result
    integer::n,mask,i,k,total
    type(linear_constraints)::cs
    type(optimizer_result)::tmp
    real(dp)::best
    real(dp), allocatable :: projected(:)
    n=size(mu);best=huge(1.0_dp);call set_failure(result,n,"no feasible subset",1)
    if(n>22)then
      cs=con
      call minvariance_portfolio(mu,sigma,cs,tmp)
      if(tmp%converged)then
        do i=1,n
          if(tmp%weights(i)<min_buyin)tmp%weights(i)=0.0_dp
        end do
        call project_feasible(tmp%weights,con,projected,tmp%status)
        tmp%weights=projected
        result=tmp;result%message="cardinality heuristic for n > 22"
      end if
      return
    end if
    total=2**n-1
    do mask=1,total
      k=popcnt(mask);if(k>max_assets)cycle
      cs=con
      do i=1,n
        if(.not.btest(mask,i-1))then
          cs%lower(i)=0.0_dp
          cs%upper(i)=0.0_dp
        else
          cs%lower(i)=max(cs%lower(i),min_buyin)
        end if
      end do
      call minvariance_portfolio(mu,sigma,cs,tmp)
      if(tmp%converged .and. tmp%risk<best)then
        best=tmp%risk
        result=tmp
      end if
    end do
    if(result%converged)result%message="exact cardinality minimum variance"
  end subroutine cardinality_minvariance_portfolio

  subroutine fill_result(result,w,obj,it,eval,conv,msg)
    type(optimizer_result),intent(out)::result
    real(dp),intent(in)::w(:),obj
    integer,intent(in)::it,eval
    logical,intent(in)::conv
    character(len=*),intent(in)::msg
    allocate(result%weights(size(w)))
    result%weights=w
    result%objective=obj
    result%iterations=it
    result%evaluations=eval
    result%converged=conv;result%status=merge(0,1,conv);result%message=msg
  end subroutine fill_result

  subroutine set_failure(result,n,msg,status)
    type(optimizer_result),intent(out)::result
    integer,intent(in)::n,status
    character(len=*),intent(in)::msg
    allocate(result%weights(n))
    result%weights=0.0_dp
    result%status=status
    result%converged=.false.
    result%message=msg
  end subroutine set_failure
end module fportfolio_optimization
