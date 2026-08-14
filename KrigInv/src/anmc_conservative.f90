! SPDX-License-Identifier: GPL-3.0-only
module anmc_conservative
  use anmc_kinds, only : dp
  use anmc_types, only : conservative_result, probability_estimate, simulation_control
  use anmc_math, only : pmvnorm, probability_control, probability_result, genz_bretz
  use anmc_utils, only : normal_cdf_local, sort_indices_descending, gather_vector, gather_matrix, &
                         positive_infinity, negative_infinity
  use anmc_probabilities, only : proba_max, proba_min
  implicit none
  private
  public :: conservative_estimate

contains

  function conservative_estimate(alpha,mean,covariance,design,threshold,pn,excursion_type,algo, &
                                 sim_control,prob_control,verb) result(res)
    real(dp),intent(in),optional::alpha
    real(dp),intent(in)::mean(:),covariance(:,:),design(:,:),threshold
    real(dp),intent(in),optional::pn(:)
    character(len=*),intent(in),optional::excursion_type,algo
    type(simulation_control),intent(in),optional::sim_control
    type(probability_control),intent(in),optional::prob_control
    integer,intent(in),optional::verb
    type(conservative_result)::res

    type(probability_control)::pctl
    real(dp)::a,prod_p,proba_left,proba_right,var_left,var_right,proba_eval,var_eval,time_ganmc
    real(dp),allocatable::coverage(:),sorted_p(:),mu_s(:),sigma_s(:,:),e_s(:,:)
    integer,allocatable::ord(:),cols(:)
    integer::n,i,ind_max,left_idx,right_idx,next_eval,meth,v,qlo,qhi
    logical::is_above
    character(len=16)::alg

    a=0.95_dp; if(present(alpha)) a=alpha
    pctl=genz_bretz(); if(present(prob_control)) pctl=prob_control
    v=0; if(present(verb)) v=verb
    is_above=.true.; if(present(excursion_type)) is_above=(excursion_type(1:1)=='>')
    alg='GANMC'; if(present(algo)) alg=adjustl(algo)
    n=size(mean)
    if(n<1 .or. size(covariance,1)/=n .or. size(covariance,2)/=n .or. size(design,1)/=n) then
      res%ok=.false.; res%message='invalid mean/covariance/design dimensions'; allocate(res%set(0)); return
    end if
    allocate(coverage(n))
    if(present(pn)) then
      if(size(pn)/=n) then; res%ok=.false.; res%message='pn length mismatch'; allocate(res%set(0)); return; end if
      coverage=pn
    else
      do i=1,n
        if(covariance(i,i)>0.0_dp) then
          if(is_above) then
            coverage(i)=normal_cdf_local((mean(i)-threshold)/sqrt(covariance(i,i)))
          else
            coverage(i)=normal_cdf_local((threshold-mean(i))/sqrt(covariance(i,i)))
          end if
        else
          if(is_above) then; coverage(i)=merge(1.0_dp,0.0_dp,mean(i)>threshold)
          else; coverage(i)=merge(1.0_dp,0.0_dp,mean(i)<threshold); end if
        end if
      end do
    end if
    ord=sort_indices_descending(coverage); allocate(sorted_p(n))
    do i=1,n; sorted_p(i)=coverage(ord(i)); end do

    ! Upstream used a fixed productPn(10000) and indexed past n when all
    ! sorted probabilities exceeded alpha.  This bounded version preserves
    ! the intended stopping rule without the overrun.
    ind_max=1; prod_p=sorted_p(1); left_idx=1
    do while(ind_max<n .and. sorted_p(ind_max)>a)
      if(prod_p>a) left_idx=ind_max
      ind_max=ind_max+1
      prod_p=prod_p*sorted_p(ind_max)
    end do
    if(prod_p>a) left_idx=ind_max
    right_idx=ind_max
    left_idx=max(1,min(left_idx,right_idx))

    mu_s=gather_vector(mean,ord(1:ind_max))
    sigma_s=gather_matrix(covariance,ord(1:ind_max),ord(1:ind_max))
    cols=[(i,i=1,size(design,2))]
    e_s=gather_matrix(design,ord(1:ind_max),cols)
    meth=4

    if(ind_max<500) then
      call direct_prefix_probability(right_idx,is_above,threshold,mu_s,sigma_s,pctl,proba_right,var_right)
      call direct_prefix_probability(left_idx,is_above,threshold,mu_s,sigma_s,pctl,proba_left,var_left)
    else
      qlo=min(10,left_idx); qhi=min(20,left_idx)
      call approximate_prefix(left_idx,is_above,10.0_dp,qlo,qhi,threshold,mu_s,sigma_s,e_s,sorted_p, &
                              meth,alg,sim_control,pctl,v,proba_left,var_left)
      call approximate_prefix(right_idx,is_above,10.0_dp,qlo,qhi,threshold,mu_s,sigma_s,e_s,sorted_p, &
                              meth,alg,sim_control,pctl,v,proba_right,var_right)
    end if

    do while(proba_right<a .and. right_idx-left_idx>=2)
      next_eval=ceiling(real(right_idx+left_idx,dp)/2.0_dp)
      if(right_idx<500) then
        call direct_prefix_probability(next_eval,is_above,threshold,mu_s,sigma_s,pctl,proba_eval,var_eval)
      else
        qlo=min(15,max(1,left_idx-1)); qhi=min(25,max(1,left_idx-1)); time_ganmc=10.0_dp
        if(right_idx-left_idx<=7) then
          time_ganmc=25.0_dp
          qlo=min(20,max(1,left_idx-1)); qhi=min(30,max(1,left_idx-1))
          if(right_idx-left_idx<=5) qhi=min(40,max(1,left_idx-1))
          if(right_idx-left_idx<=3) then
            time_ganmc=30.0_dp
            qlo=min(25,max(1,left_idx-1)); qhi=min(60,max(1,left_idx-1))
          end if
        end if
        call approximate_prefix(next_eval,is_above,time_ganmc,qlo,qhi,threshold,mu_s,sigma_s,e_s,sorted_p, &
                                meth,alg,sim_control,pctl,v,proba_eval,var_eval)
      end if
      if(proba_eval>a) then
        left_idx=next_eval; proba_left=proba_eval; var_left=var_eval
      else
        right_idx=next_eval; proba_right=proba_eval; var_right=var_eval
      end if
    end do

    allocate(res%set(n)); res%set=coverage>sorted_p(right_idx)
    res%level=sorted_p(right_idx); res%probability=proba_left; res%uncertainty=var_left
    res%ok=.true.; res%message='ok'
    if(v>0) write(*,'(a,f10.6,a,i0)') 'Conservative level = ',res%level,', retained points = ',count(res%set)
  end function conservative_estimate

  subroutine direct_prefix_probability(k,is_above,threshold,mu,sigma,pctl,p,var_like)
    integer,intent(in)::k
    logical,intent(in)::is_above
    real(dp),intent(in)::threshold,mu(:),sigma(:,:)
    type(probability_control),intent(in)::pctl
    real(dp),intent(out)::p,var_like
    type(probability_result)::pr
    real(dp),allocatable::lo(:),up(:)
    integer,allocatable::idx(:)
    integer::i
    allocate(lo(k),up(k)); idx=[(i,i=1,k)]
    if(is_above) then; lo=threshold; up=positive_infinity()
    else; lo=negative_infinity(); up=threshold; end if
    pr=pmvnorm(lo,up,mu(1:k),gather_matrix(sigma,idx,idx),pctl)
    p=pr%value
    ! Preserve upstream conservativeEstimate semantics: this branch stores
    ! mvtnorm's absolute error in the field called `vars`.
    var_like=pr%error
  end subroutine direct_prefix_probability

  subroutine approximate_prefix(k,is_above,budget,qlo,qhi,threshold,mu,sigma,e,sorted_p,meth,alg,sctl,pctl,v,p,var_like)
    integer,intent(in)::k,qlo,qhi,meth,v
    logical,intent(in)::is_above
    real(dp),intent(in)::budget,threshold,mu(:),sigma(:,:),e(:,:),sorted_p(:)
    character(len=*),intent(in)::alg
    type(simulation_control),intent(in),optional::sctl
    type(probability_control),intent(in)::pctl
    real(dp),intent(out)::p,var_like
    type(probability_estimate)::pe
    integer,allocatable::idx(:),cols(:)
    real(dp),allocatable::pk(:),sk(:,:),ek(:,:),pn_k(:)
    integer::i
    idx=[(i,i=1,k)]; cols=[(i,i=1,size(e,2))]
    pk=mu(1:k); sk=gather_matrix(sigma,idx,idx); ek=gather_matrix(e,idx,cols)
    if(is_above) then
      pe=proba_min(budget,threshold,pk,sk,ek,q_limits=[qlo,qhi],method=meth,algo=alg, &
                   sim_control=sctl,prob_control=pctl,verb=max(0,v-1))
    else
      pn_k=sorted_p(1:k)
      pe=proba_max(budget,threshold,pk,sk,ek,q_limits=[qlo,qhi],pn=pn_k,method=meth, &
                   algo=alg,sim_control=sctl,prob_control=pctl,verb=max(0,v-1))
    end if
    p=1.0_dp-pe%probability; var_like=pe%variance
  end subroutine approximate_prefix

end module anmc_conservative
