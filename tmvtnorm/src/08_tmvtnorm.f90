! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm
  use mvtnorm_kinds, only : dp
  use tmvtnorm_distributions, only : ptmvnorm, dtmvnorm_one, dtmvnorm, ptmvt, dtmvt_one, dtmvt
  use tmvtnorm_marginals, only : dtmvnorm_marginal, dtmvnorm_marginal2, ptmvnorm_marginal, &
    qtmvnorm_marginal, ptmvt_marginal, tail_lower, tail_upper, tail_both
  use tmvtnorm_moments, only : tmvnorm_moments_t, mtmvnorm
  use tmvtnorm_random, only : rtnorm, rtmvnorm_rejection, rtmvnorm_gibbs, rtmvnorm_gibbs_precision, &
    rtmvnorm_gibbs_linear, rtmvnorm_sparse_csc, rtmvnorm_sparse_triplet, rtmvt_rejection, rtmvt_gibbs
  use tmvtnorm_estimation, only : tmvnorm_fit_t, tmvnorm_gmm_fit_t, mle_tmvnorm, gmm_tmvnorm, &
    gmm_moments_manjunath_wilhelm, gmm_moments_lee, gmm_mw, gmm_lee
  implicit none
  private
  public :: dp
  public :: ptmvnorm,dtmvnorm_one,dtmvnorm,ptmvt,dtmvt_one,dtmvt
  public :: dtmvnorm_marginal,dtmvnorm_marginal2,ptmvnorm_marginal,qtmvnorm_marginal,ptmvt_marginal
  public :: tail_lower,tail_upper,tail_both,tmvnorm_moments_t,mtmvnorm
  public :: rtnorm,rtmvnorm_rejection,rtmvnorm_gibbs,rtmvnorm_gibbs_precision,rtmvnorm_gibbs_linear
  public :: rtmvnorm_sparse_csc,rtmvnorm_sparse_triplet,rtmvt_rejection,rtmvt_gibbs
  public :: tmvnorm_fit_t,tmvnorm_gmm_fit_t,mle_tmvnorm,gmm_tmvnorm
  public :: gmm_moments_manjunath_wilhelm,gmm_moments_lee,gmm_mw,gmm_lee
  public :: rtmvnorm,rtmvnorm2,rtmvt,algorithm_rejection,algorithm_gibbs
  integer,parameter :: algorithm_rejection=1,algorithm_gibbs=2

contains

  function rtmvnorm(n,mean,sigma,lower,upper,algorithm,dmat,burnin,thinning,start,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in),optional::algorithm,burnin,thinning,seed
    real(dp),intent(in),optional::dmat(:,:),start(:)
    real(dp),allocatable::x(:,:)
    integer::alg
    alg=algorithm_gibbs
    if(present(algorithm)) alg=algorithm
    if(present(dmat)) then
      if(alg==algorithm_rejection) then
        x=rtmvnorm_rejection(n,mean,sigma,lower,upper,dmat,seed)
      else
        x=rtmvnorm_gibbs_linear(n,mean,sigma,dmat,lower,upper,burnin,thinning,start,seed)
      end if
    else
      if(alg==algorithm_rejection) then
        x=rtmvnorm_rejection(n,mean,sigma,lower,upper,seed=seed)
      else
        x=rtmvnorm_gibbs(n,mean,sigma,lower,upper,burnin,thinning,start,seed)
      end if
    end if
  end function rtmvnorm

  function rtmvnorm2(n,mean,sigma,dmat,lower,upper,algorithm,burnin,thinning,start,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),dmat(:,:),lower(:),upper(:)
    integer,intent(in),optional::algorithm,burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::x(:,:)
    integer::alg
    alg=algorithm_gibbs
    if(present(algorithm)) alg=algorithm
    if(alg==algorithm_rejection) then
      x=rtmvnorm_rejection(n,mean,sigma,lower,upper,dmat,seed)
    else
      x=rtmvnorm_gibbs_linear(n,mean,sigma,dmat,lower,upper,burnin,thinning,start,seed)
    end if
  end function rtmvnorm2

  function rtmvt(n,mean,sigma,df,lower,upper,algorithm,burnin,thinning,start,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),df,lower(:),upper(:)
    integer,intent(in),optional::algorithm,burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::x(:,:)
    integer::alg
    alg=algorithm_rejection
    if(present(algorithm)) alg=algorithm
    if(alg==algorithm_rejection) then
      x=rtmvt_rejection(n,mean,sigma,df,lower,upper,seed)
    else
      x=rtmvt_gibbs(n,mean,sigma,df,lower,upper,burnin,thinning,start,seed)
    end if
  end function rtmvt

end module tmvtnorm
