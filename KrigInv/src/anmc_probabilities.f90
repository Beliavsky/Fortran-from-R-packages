! SPDX-License-Identifier: GPL-3.0-only
module anmc_probabilities
  use anmc_kinds, only : dp
  use anmc_types, only : probability_estimate, anmc_problem, simulation_control, active_dims_result
  use anmc_math, only : pmvnorm, probability_control, probability_result, genz_bretz, cholesky_lower, inverse_spd
  use anmc_utils, only : normal_cdf_local, complement_indices, gather_vector, gather_matrix, &
                         positive_infinity, negative_infinity, seed_fortran_rng
  use anmc_active, only : select_active_dims, select_q_dims
  use anmc_mc, only : mc_gauss, anmc_gauss
  implicit none
  private
  public :: proba_max, proba_min

contains

  function proba_max(c_bdg,threshold,mu,sigma,e,q,q_limits,pn,method,algo,sim_control,prob_control,seed,verb) result(res)
    real(dp),intent(in)::c_bdg,threshold,mu(:),sigma(:,:)
    real(dp),intent(in),optional::e(:,:),pn(:)
    integer,intent(in),optional::q,q_limits(2),method,seed,verb
    character(len=*),intent(in),optional::algo
    type(simulation_control),intent(in),optional::sim_control
    type(probability_control),intent(in),optional::prob_control
    type(probability_estimate)::res
    res=proba_core(.true.,c_bdg,threshold,mu,sigma,e,q,q_limits,pn,method,algo,sim_control,prob_control,seed,verb)
  end function proba_max

  function proba_min(c_bdg,threshold,mu,sigma,e,q,q_limits,pn,method,algo,sim_control,prob_control,seed,verb) result(res)
    real(dp),intent(in)::c_bdg,threshold,mu(:),sigma(:,:)
    real(dp),intent(in),optional::e(:,:),pn(:)
    integer,intent(in),optional::q,q_limits(2),method,seed,verb
    character(len=*),intent(in),optional::algo
    type(simulation_control),intent(in),optional::sim_control
    type(probability_control),intent(in),optional::prob_control
    type(probability_estimate)::res
    res=proba_core(.false.,c_bdg,threshold,mu,sigma,e,q,q_limits,pn,method,algo,sim_control,prob_control,seed,verb)
  end function proba_min

  function proba_core(is_max,c_bdg,threshold,mu,sigma,e,q,q_limits,pn,method,algo,sim_control,prob_control,seed,verb) result(res)
    logical,intent(in)::is_max
    real(dp),intent(in)::c_bdg,threshold,mu(:),sigma(:,:)
    real(dp),intent(in),optional::e(:,:),pn(:)
    integer,intent(in),optional::q,q_limits(2),method,seed,verb
    character(len=*),intent(in),optional::algo
    type(simulation_control),intent(in),optional::sim_control
    type(probability_control),intent(in),optional::prob_control
    type(probability_estimate)::res

    type(probability_control)::pctl
    type(probability_result)::pr
    type(anmc_problem)::problem
    type(active_dims_result)::asel
    real(dp),allocatable::design(:,:),coverage(:),weights(:),mu_eq(:),k_eq(:,:),lower(:),upper(:)
    real(dp),allocatable::l(:,:),invxx(:,:),sigma_xy(:,:),sigma_yy(:,:),sigma_cond(:,:)
    integer,allocatable::ind_q(:),ind_c(:),cols(:)
    integer::n,meth,attempt,v
    logical::ok,use_anmc
    character(len=256)::msg
    character(len=16)::alg
    real(dp)::pprime,var_pprime

    n=size(mu); meth=4; if(present(method)) meth=method
    v=0; if(present(verb)) v=verb
    pctl=genz_bretz(); if(present(prob_control)) pctl=prob_control
    if(present(seed)) call seed_fortran_rng(seed)
    if(size(sigma,1)/=n .or. size(sigma,2)/=n .or. n<1) then
      res%ok=.false.; res%message='invalid mean/covariance dimensions'; return
    end if
    if(present(e)) then
      if(size(e,1)/=n) then; res%ok=.false.; res%message='design row count must equal length(mu)'; return; end if
      design=e
    else
      allocate(design(n,1))
      if(n==1) then; design(1,1)=0.0_dp
      else
        design(:,1)=[(real(attempt-1,dp)/real(n-1,dp),attempt=1,n)]
      end if
    end if
    allocate(coverage(n),weights(n))
    if(present(pn)) then
      if(size(pn)/=n) then; res%ok=.false.; res%message='pn length mismatch'; return; end if
      coverage=pn
    else
      do attempt=1,n
        if(sigma(attempt,attempt)>0.0_dp) then
          if(is_max) then
            coverage(attempt)=normal_cdf_local((threshold-mu(attempt))/sqrt(sigma(attempt,attempt)))
          else
            coverage(attempt)=normal_cdf_local((mu(attempt)-threshold)/sqrt(sigma(attempt,attempt)))
          end if
        else
          if(is_max) then
            coverage(attempt)=merge(1.0_dp,0.0_dp,mu(attempt)<=threshold)
          else
            coverage(attempt)=merge(1.0_dp,0.0_dp,mu(attempt)>=threshold)
          end if
        end if
      end do
    end if
    weights=1.0_dp-coverage

    if(present(q)) then
      ind_q=select_active_dims(max(1,min(q,n)),design,threshold,mu,sigma,weights,meth)
    else if(present(q_limits)) then
      asel=select_q_dims(design,threshold,mu,sigma,weights,meth,q_limits,pctl,.true.)
      ind_q=asel%ind_q
    else
      asel=select_q_dims(design,threshold,mu,sigma,weights,meth,prob_control=pctl,reduced_return=.true.)
      ind_q=asel%ind_q
    end if
    if(size(ind_q)==0) then; res%ok=.false.; res%message='active-dimension selection failed'; return; end if

    ! Re-select if the active covariance is singular, as the R code does.
    do attempt=1,100
      mu_eq=gather_vector(mu,ind_q); k_eq=gather_matrix(sigma,ind_q,ind_q)
      call cholesky_lower(k_eq,l,ok,msg)
      if(ok) exit
      if(present(q)) then
        ind_q=select_active_dims(size(ind_q),design,threshold,mu,sigma,weights,meth)
      else
        ind_q=select_active_dims(size(ind_q),design,threshold,mu,sigma,weights,meth)
      end if
    end do
    if(.not.ok) then; res%ok=.false.; res%message='active covariance remained non-positive-definite'; return; end if

    allocate(lower(size(ind_q)),upper(size(ind_q)))
    if(is_max) then
      lower=negative_infinity(); upper=threshold
    else
      lower=threshold; upper=positive_infinity()
    end if
    pr=pmvnorm(lower,upper,mu_eq,k_eq,pctl)
    pprime=1.0_dp-pr%value
    res%pq=pprime; res%pq_error=pr%error; res%q=size(ind_q); res%ind_q=ind_q
    cols=[(attempt,attempt=1,size(design,2))]
    res%eq=gather_matrix(design,ind_q,cols)
    if((1.0_dp-pprime)<pr%error .or. size(ind_q)==n) then
      res%probability=pprime
      res%variance=(2.0_dp/7.0_dp*pr%error)**2
      res%rq=0.0_dp; res%has_remainder=.false.
      res%ok=(pr%inform==0 .or. pr%inform==1); res%message=trim(pr%message)
      return
    end if

    ind_c=complement_indices(n,ind_q)
    problem%mu_eq=mu_eq; problem%sigma_eq=k_eq; problem%threshold=threshold
    problem%mu_emq=gather_vector(mu,ind_c)
    call inverse_spd(k_eq,invxx,ok,msg)
    if(.not.ok) then; res%ok=.false.; res%message=msg; return; end if
    sigma_xy=gather_matrix(sigma,ind_q,ind_c)
    problem%ww_cond_q=matmul(transpose(sigma_xy),invxx)
    sigma_yy=gather_matrix(sigma,ind_c,ind_c)
    sigma_cond=sigma_yy-matmul(problem%ww_cond_q,sigma_xy)
    call cholesky_lower(sigma_cond,l,ok,msg)
    if(.not.ok) then; res%ok=.false.; res%message='conditional covariance Cholesky failed: '//trim(msg); return; end if
    problem%sigma_cond_q_chol=transpose(l)

    alg='ANMC'; if(present(algo)) alg=adjustl(algo)
    use_anmc=(trim(alg)=='ANMC' .or. trim(alg)=='GANMC')
    if(use_anmc) then
      if(is_max) then
        res%remainder=anmc_gauss(c_bdg,problem,0.45_dp,'M',sim_control,pctl,v)
      else
        res%remainder=anmc_gauss(c_bdg,problem,0.45_dp,'m',sim_control,pctl,v)
      end if
    else
      if(is_max) then
        res%remainder=mc_gauss(c_bdg,problem,0.20_dp,'M',sim_control=sim_control,prob_control=pctl,verb=v)
      else
        res%remainder=mc_gauss(c_bdg,problem,0.20_dp,'m',sim_control=sim_control,prob_control=pctl,verb=v)
      end if
    end if
    if(.not.res%remainder%ok) then; res%ok=.false.; res%message=res%remainder%message; return; end if
    res%rq=res%remainder%estim
    res%probability=pprime+res%rq*(1.0_dp-pprime)
    var_pprime=(pr%error/3.0_dp)**2
    res%variance=(1.0_dp-res%rq)**2*var_pprime + res%remainder%var_est*(1.0_dp-pprime)**2 + &
                 res%remainder%var_est*var_pprime
    res%has_remainder=.true.; res%ok=.true.; res%message='ok'
  end function proba_core

end module anmc_probabilities
