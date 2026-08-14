! SPDX-License-Identifier: GPL-3.0-only
module anmc_mc
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use anmc_kinds, only : dp
  use anmc_types, only : anmc_problem, mc_result, mc_params, simulation_control
  use anmc_math, only : probability_control, genz_bretz
  use anmc_utils, only : wall_time_seconds, mean_value, sample_variance, mean_finite, &
                         linear_fit, slope_through_origin, positive_infinity, negative_infinity
  use anmc_sampling, only : mvrnorm_arma, trmvrnorm_rej_cpp
  implicit none
  private
  public :: mc_gauss, anmc_gauss

contains

  function mc_gauss(comp_bdg, problem, delta, excursion_type, params, sim_control, prob_control, verb) result(res)
    real(dp), intent(in) :: comp_bdg
    type(anmc_problem), intent(in) :: problem
    real(dp), intent(in), optional :: delta
    character(len=*), intent(in), optional :: excursion_type
    type(mc_params), intent(in), optional :: params
    type(simulation_control), intent(in), optional :: sim_control
    type(probability_control), intent(in), optional :: prob_control
    integer, intent(in), optional :: verb
    type(mc_result) :: res

    type(simulation_control) :: sctl
    type(probability_control) :: pctl
    real(dp) :: del, time_start, time1x, time1y, t0, cx0, alpha_cost, beta0, t_eval_g
    real(dp), allocatable :: upper(:), lower(:), sims_x(:,:), temp(:,:), mu_y(:), y(:,:), g(:)
    real(dp), allocatable :: tt_x(:), ii_x(:), tt_y(:), ii_y(:)
    integer :: size_x, size_y, i, j, n0, nstar, step_x, step_y, n_batch, v, n_done
    logical :: is_max, ok
    character(len=256) :: msg

    sctl = simulation_control(); if (present(sim_control)) sctl = sim_control
    pctl = genz_bretz(); if (present(prob_control)) pctl = prob_control
    del = 0.1_dp; if (present(delta)) del = delta
    v = 0; if (present(verb)) v = verb
    is_max = .true.
    if (present(excursion_type)) is_max = (excursion_type(1:1) /= 'm')
    size_x = size(problem%mu_eq); size_y = size(problem%mu_emq)
    allocate(upper(size_x), lower(size_x))
    if (is_max) then
      upper = problem%threshold; lower = negative_infinity()
    else
      upper = positive_infinity(); lower = problem%threshold
    end if

    time_start = wall_time_seconds()
    if (.not. present(params)) then
      ! Estimate the per-outer-sample cost, matching lm(ttX ~ ii + 0).
      t0 = wall_time_seconds()
      sims_x = trmvrnorm_rej_cpp(1,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
      time1x = max(wall_time_seconds()-t0, 1.0e-9_dp)
      if (.not. ok) then
        res%ok=.false.; res%message='truncated-normal calibration failed'; return
      end if
      allocate(tt_x(20),ii_x(20))
      step_x = max(1, floor((comp_bdg*del*0.4_dp/time1x - 20.0_dp)/190.0_dp))
      step_x = min(step_x, max(1,sctl%max_outer/210))
      do i=1,20
        n_batch = 1 + (i-1)*step_x
        ii_x(i)=real(n_batch,dp)
        t0=wall_time_seconds()
        temp=trmvrnorm_rej_cpp(n_batch,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
        tt_x(i)=max(wall_time_seconds()-t0,0.0_dp)*1.03_dp
        if(.not.ok) then; res%ok=.false.; res%message='truncated-normal calibration failed'; return; end if
        call append_columns(sims_x,temp)
      end do
      cx0=max(slope_through_origin(ii_x,tt_x),1.0e-12_dp)

      allocate(mu_y(size_y))
      t0=wall_time_seconds()
      mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,1)-problem%mu_eq)
      y=mvrnorm_arma(1,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
      time1y=max(wall_time_seconds()-t0,1.0e-9_dp)
      if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
      allocate(tt_y(20),ii_y(20))
      step_y=max(1,floor((comp_bdg*del*0.5_dp/time1y-20.0_dp)/190.0_dp))
      step_y=min(step_y,max(1,sctl%max_inner/210))
      do i=1,20
        n_batch=1+(i-1)*step_y
        ii_y(i)=real(n_batch,dp)
        t0=wall_time_seconds()
        mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,1)-problem%mu_eq)
        y=mvrnorm_arma(n_batch,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
        tt_y(i)=max(wall_time_seconds()-t0,0.0_dp)*1.03_dp
        if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
      end do
      call linear_fit(ii_y(2:20),tt_y(2:20),alpha_cost,beta0)
      alpha_cost=max(alpha_cost,0.0_dp); beta0=max(beta0,1.0e-12_dp)
      t_eval_g=estimate_g_cost(is_max,problem%threshold,size_y)
      n0=size(sims_x,2)
    else
      cx0=max(params%cx,1.0e-12_dp)
      alpha_cost=max(params%alpha_cost,0.0_dp)
      beta0=max(params%beta,1.0e-12_dp)
      t_eval_g=max(params%eval_g,0.0_dp)
      n0=0
      allocate(sims_x(size_x,0))
    end if

    nstar=nint(comp_bdg/max(cx0+beta0+t_eval_g,1.0e-12_dp))
    nstar=max(1,min(nstar,sctl%max_outer))
    if (.not. present(params)) nstar=max(nstar,n0)
    if (present(params)) then
      sims_x=trmvrnorm_rej_cpp(nstar,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
      if(.not.ok) then; res%ok=.false.; res%message='truncated-normal production sampling failed'; return; end if
    else if(nstar>n0) then
      temp=trmvrnorm_rej_cpp(nstar-n0,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
      if(.not.ok) then; res%ok=.false.; res%message='truncated-normal production sampling failed'; return; end if
      call append_columns(sims_x,temp)
    else
      nstar=n0
    end if

    allocate(g(nstar))
    if (.not. allocated(mu_y)) allocate(mu_y(size_y))
    g=0.0_dp
    n_done=0
    do j=1,nstar
      mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,j)-problem%mu_eq)
      y=mvrnorm_arma(1,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
      if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
      g(j)=merge(1.0_dp,0.0_dp,excursion(y(:,1),problem%threshold,is_max))
      n_done=j
      if (sctl%enforce_budget .and. modulo(j,100)==0) then
        if (wall_time_seconds()-time_start >= comp_bdg*sctl%time_guard) exit
      end if
    end do
    nstar=n_done
    if(nstar<size(g)) g=g(1:nstar)
    res%estim=mean_value(g)
    if(nstar>1) res%var_est=sample_variance(g)/real(nstar,dp)
    res%params%n=nstar; res%params%m=1; res%params%cx=cx0; res%params%alpha_cost=alpha_cost
    res%params%beta=beta0; res%params%eval_g=t_eval_g
    res%time_part1=max(0.0_dp,wall_time_seconds()-time_start)
    res%time_total=res%time_part1
    if(v>0) write(*,'(a,f10.6,a,i0)') 'MC estimate = ',res%estim,', n = ',nstar
  end function mc_gauss

  function anmc_gauss(comp_bdg, problem, delta, excursion_type, sim_control, prob_control, verb, &
                      fixed_n, fixed_m) result(res)
    real(dp), intent(in) :: comp_bdg
    type(anmc_problem), intent(in) :: problem
    real(dp), intent(in), optional :: delta
    character(len=*), intent(in), optional :: excursion_type
    type(simulation_control), intent(in), optional :: sim_control
    type(probability_control), intent(in), optional :: prob_control
    integer, intent(in), optional :: verb, fixed_n, fixed_m
    type(mc_result) :: res

    type(simulation_control) :: sctl
    type(probability_control) :: pctl
    real(dp) :: del,time_start,time1x,time1y,t0,cx0,cx,alpha_cost,beta0,t_eval_g,c_adj
    real(dp) :: ratio,mstar_real,eps_m,time_part1,var_g,var_g0,var_exp
    real(dp), allocatable :: upper(:),lower(:),sims_x(:,:),temp(:,:),mu_y(:),y(:,:)
    real(dp), allocatable :: ii(:),tt(:),g0(:,:),exp0(:),var0(:),ratios(:),gcur(:),expfull(:),varfull(:)
    integer :: size_x,size_y,n_tests,step_i,n0,m0,n_betas,nstar,mstar,i,j,indm,v,n_batch
    logical :: is_max,ok
    character(len=256) :: msg

    sctl=simulation_control(); if(present(sim_control)) sctl=sim_control
    pctl=genz_bretz(); if(present(prob_control)) pctl=prob_control
    del=0.4_dp; if(present(delta)) del=delta
    v=0; if(present(verb)) v=verb
    is_max=.true.; if(present(excursion_type)) is_max=(excursion_type(1:1)/='m')
    size_x=size(problem%mu_eq); size_y=size(problem%mu_emq)
    allocate(upper(size_x),lower(size_x))
    if(is_max) then; upper=problem%threshold; lower=negative_infinity()
    else; upper=positive_infinity(); lower=problem%threshold; end if
    time_start=wall_time_seconds()

    t0=wall_time_seconds()
    sims_x=trmvrnorm_rej_cpp(1,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
    time1x=max(wall_time_seconds()-t0,1.0e-9_dp)
    if(.not.ok) then; res%ok=.false.; res%message='truncated-normal calibration failed'; return; end if
    n_tests=max(2,min(floor(comp_bdg*del*0.4_dp/time1x),floor(2060.0_dp/61.0_dp)))
    n_tests=min(n_tests,33)
    step_i=max(1,ceiling(120.0_dp/real(n_tests,dp)))
    n_tests=1+(120-1)/step_i
    allocate(ii(n_tests),tt(n_tests))
    do i=1,n_tests
      n_batch=1+(i-1)*step_i
      ii(i)=real(n_batch,dp)
      t0=wall_time_seconds()
      temp=trmvrnorm_rej_cpp(n_batch,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
      tt(i)=max(wall_time_seconds()-t0,0.0_dp)*1.03_dp
      if(.not.ok) then; res%ok=.false.; res%message='truncated-normal calibration failed'; return; end if
      call append_columns(sims_x,temp)
    end do
    call linear_fit(ii,tt,cx0,cx)
    cx=max(cx,1.0e-12_dp); cx0=max(cx0,0.0_dp)

    allocate(mu_y(size_y))
    t0=wall_time_seconds()
    mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,1)-problem%mu_eq)
    alpha_cost=max(wall_time_seconds()-t0,0.0_dp)
    y=mvrnorm_arma(1,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
    time1y=max(wall_time_seconds()-t0,1.0e-9_dp)
    if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
    n_betas=max(1,ceiling(comp_bdg*del*0.1_dp/max(time1y-alpha_cost,1.0e-9_dp)))
    n_betas=min(n_betas,sctl%max_inner)
    t0=wall_time_seconds(); y=mvrnorm_arma(n_betas,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
    beta0=max((wall_time_seconds()-t0)/real(n_betas,dp),1.0e-12_dp)
    if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
    t_eval_g=estimate_g_cost(is_max,problem%threshold,size_y)
    c_adj=comp_bdg*del-time1x-sum(tt)-alpha_cost-beta0*real(n_betas,dp)-time1y-154.0_dp*t_eval_g
    n0=size(sims_x,2)
    m0=max(30,floor(((c_adj-cx0)/real(max(1,n0),dp)-cx-alpha_cost)/max(beta0+t_eval_g,1.0e-12_dp)))
    m0=min(m0,sctl%max_inner)
    if(present(fixed_m)) m0=max(1,min(fixed_m,sctl%max_inner))

    allocate(g0(m0,n0),exp0(n0),var0(n0),gcur(m0))
    do j=1,n0
      mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,j)-problem%mu_eq)
      y=mvrnorm_arma(m0,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
      if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
      do i=1,m0
        g0(i,j)=merge(1.0_dp,0.0_dp,excursion(y(:,i),problem%threshold,is_max))
      end do
      exp0(j)=mean_value(g0(:,j)); var0(j)=sample_variance(g0(:,j))
    end do
    res%estim0=mean_value(exp0)
    var_exp=sample_variance(exp0)
    allocate(ratios(m0)); ratios=0.0_dp
    if(var_exp>tiny(1.0_dp)) then
      do i=1,m0
        ratios(i)=sample_variance(g0(i,:))/var_exp-1.0_dp
      end do
      ratio=mean_finite(ratios,0.25_dp)
    else
      ratio=0.25_dp
    end if
    if(.not.ieee_is_finite(ratio)) ratio=0.25_dp
    ratio=max(ratio,0.0_dp)
    mstar_real=sqrt(max(0.0_dp,((cx0+cx+3.0_dp*alpha_cost)*ratio)/max(beta0+t_eval_g,1.0e-12_dp)))
    if(.not.ieee_is_finite(mstar_real)) mstar_real=real(m0,dp)
    eps_m=mstar_real-floor(mstar_real)
    if(eps_m<0.5_dp*(2.0_dp*mstar_real+1.0_dp-sqrt(4.0_dp*mstar_real*mstar_real+1.0_dp))) then
      mstar=max(floor(mstar_real),1)
    else
      mstar=max(ceiling(mstar_real),1)
    end if
    mstar=min(mstar,sctl%max_inner)
    if(present(fixed_m)) mstar=max(1,min(fixed_m,sctl%max_inner))
    time_part1=wall_time_seconds()-time_start
    nstar=nint((comp_bdg-0.9_dp*time_part1-cx0)/max(cx+alpha_cost+beta0*real(mstar,dp),1.0e-12_dp))
    nstar=max(n0,min(nstar,sctl%max_outer))
    if(present(fixed_n)) nstar=max(n0,min(fixed_n,sctl%max_outer))
    if(nstar>n0) then
      temp=trmvrnorm_rej_cpp(nstar-n0,problem%mu_eq,problem%sigma_eq,lower,upper,v-1,sctl,pctl,ok)
      if(.not.ok) then; res%ok=.false.; res%message='truncated-normal production sampling failed'; return; end if
      call append_columns(sims_x,temp)
    end if

    deallocate(gcur); allocate(gcur(mstar),expfull(nstar),varfull(nstar))
    indm=min(m0,mstar)
    do j=1,nstar
      gcur=0.0_dp
      if(j<=n0) then
        gcur(1:indm)=g0(1:indm,j)
        if(mstar>m0) then
          ! Upstream tested indM > m0, which can never be true because
          ! indM=min(m0,mStar).  The intended condition is mStar > m0.
          mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,j)-problem%mu_eq)
          y=mvrnorm_arma(mstar-m0,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
          if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
          do i=m0+1,mstar
            gcur(i)=merge(1.0_dp,0.0_dp,excursion(y(:,i-m0),problem%threshold,is_max))
          end do
        end if
      else
        mu_y=problem%mu_emq+matmul(problem%ww_cond_q,sims_x(:,j)-problem%mu_eq)
        y=mvrnorm_arma(mstar,mu_y,problem%sigma_cond_q_chol,1,ok,msg)
        if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
        do i=1,mstar
          gcur(i)=merge(1.0_dp,0.0_dp,excursion(y(:,i),problem%threshold,is_max))
        end do
      end if
      expfull(j)=mean_value(gcur)
      if(mstar>1) then; varfull(j)=sample_variance(gcur); else; varfull(j)=0.0_dp; end if
    end do
    res%estim=mean_value(expfull)
    var_g=sample_variance(expfull)+mean_value(varfull)
    res%var_est=var_g/real(nstar,dp)-real(mstar-1,dp)/real(nstar*mstar,dp)*mean_value(varfull)
    var_g0=sample_variance(exp0)+mean_value(var0)
    res%var_est0=var_g0/real(n0,dp)-real(m0-1,dp)/real(n0*m0,dp)*mean_value(var0)
    res%ratio0=ratio; res%n0=n0; res%m0=m0
    res%params%n=nstar; res%params%m=mstar; res%params%cx=cx; res%params%cx0=cx0
    res%params%alpha_cost=alpha_cost; res%params%beta=beta0; res%params%eval_g=t_eval_g
    res%exp_y_cond_x=expfull; res%var_y_cond_x=varfull
    res%time_part1=time_part1; res%time_total=wall_time_seconds()-time_start
    if(v>0) write(*,'(a,f10.6,a,i0,a,i0)') 'ANMC estimate = ',res%estim,', n = ',nstar,', m = ',mstar
  end function anmc_gauss

  logical function excursion(x,threshold,is_max) result(hit)
    real(dp),intent(in)::x(:),threshold
    logical,intent(in)::is_max
    if(size(x)==0) then
      hit=.false.
    else if(is_max) then
      hit=maxval(x)>threshold
    else
      hit=minval(x)<threshold
    end if
  end function excursion

  real(dp) function estimate_g_cost(is_max,threshold,n) result(cost)
    logical,intent(in)::is_max
    real(dp),intent(in)::threshold
    integer,intent(in)::n
    real(dp)::t0,dummy(154),x(max(1,n))
    integer::i
    logical::h
    x=threshold
    do i=1,154
      t0=wall_time_seconds(); h=excursion(x,threshold,is_max)
      dummy(i)=max(wall_time_seconds()-t0,0.0_dp)
      if(h) dummy(i)=dummy(i)+0.0_dp
    end do
    call sort_real(dummy)
    ! Upstream drops extremes and uses the 0.99 quantile.  A high order
    ! statistic is a stable equivalent for such tiny timing intervals.
    cost=max(dummy(150),0.0_dp)
  end function estimate_g_cost

  subroutine append_columns(a,b)
    real(dp),allocatable,intent(inout)::a(:,:)
    real(dp),intent(in)::b(:,:)
    real(dp),allocatable::tmp(:,:)
    integer::na,nb
    if(size(a,1)/=size(b,1)) error stop 'append_columns: row mismatch'
    na=size(a,2); nb=size(b,2)
    allocate(tmp(size(a,1),na+nb))
    if(na>0) tmp(:,1:na)=a
    if(nb>0) tmp(:,na+1:na+nb)=b
    call move_alloc(tmp,a)
  end subroutine append_columns

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

end module anmc_mc
