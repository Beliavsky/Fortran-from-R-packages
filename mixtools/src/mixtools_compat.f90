! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_compat
  use mixtools_kinds, only : dp
  use mixtools_types
  use mixtools_rng, only : rng_state
  use mixtools_distributions, only : rnormmix, rmvnormmix
  use mixtools_utilities, only : ellipse_points, permutations
  use mixtools_parametric
  use mixtools_regression
  use mixtools_support, only : posterior_beta_intervals
  use mixtools_diagnostics, only : test_equality_normal, test_equality_regression
  implicit none
  private
  public :: normalmix_init, mvnormalmix_init, gammamix_init, multmix_init
  public :: regmix_init, regmix_lambda_init, repnormmix_init
  public :: logisregmix_init, poisregmix_init, segregmix_init, flaremix_init
  public :: regmix_mixed_init, try_flare, post_beta, regcr, test_equality
  public :: ellipse, perm, normmix_sim, normmixrm_sim
  interface test_equality
    module procedure test_equality_vector
    module procedure test_equality_matrix
  end interface test_equality
contains
  function one_step_control(control) result(ctl)
    type(em_control),intent(in),optional::control
    type(em_control)::ctl
    ctl=em_control();if(present(control))ctl=control;ctl%max_iterations=1
  end function one_step_control

  subroutine ellipse(mu,sigma,alpha,npoints,points,status)
    real(dp),intent(in)::mu(2),sigma(2,2),alpha
    integer,intent(in)::npoints
    real(dp),allocatable,intent(out)::points(:,:)
    integer,intent(out)::status
    call ellipse_points(mu,sigma,alpha,npoints,points,status)
  end subroutine ellipse

  subroutine perm(n,r,out,status)
    integer,intent(in)::n,r
    integer,allocatable,intent(out)::out(:,:)
    integer,intent(out)::status
    call permutations(n,r,out,status)
  end subroutine perm

  subroutine normmix_sim(rng,n,lambda,mu,sigma,x,component)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),mu(:),sigma(:)
    real(dp),intent(out)::x(n)
    integer,intent(out),optional::component(n)
    call rnormmix(rng,n,lambda,mu,sigma,x,component)
  end subroutine normmix_sim

  subroutine normmixrm_sim(rng,n,lambda,mu,sigma,x,component,status)
    type(rng_state),intent(inout)::rng
    integer,intent(in)::n
    real(dp),intent(in)::lambda(:),mu(:,:),sigma(:,:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer,allocatable,intent(out),optional::component(:)
    integer,intent(out)::status
    call rmvnormmix(rng,n,lambda,mu,sigma,x,component,status)
  end subroutine normmixrm_sim

  subroutine normalmix_init(x,k,result,control)
    real(dp),intent(in)::x(:);integer,intent(in)::k;type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call normalmix_em(x,k,result,one_step_control(control))
  end subroutine normalmix_init

  subroutine mvnormalmix_init(x,k,result,control)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k;type(mv_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call mvnormalmix_em(x,k,result,one_step_control(control))
  end subroutine mvnormalmix_init

  subroutine gammamix_init(x,k,result,control)
    real(dp),intent(in)::x(:);integer,intent(in)::k;type(gamma_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call gammamix_em(x,k,result,one_step_control(control))
  end subroutine gammamix_init

  subroutine multmix_init(y,k,result,control)
    real(dp),intent(in)::y(:,:);integer,intent(in)::k;type(multinomial_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call multmix_em(y,k,result,one_step_control(control))
  end subroutine multmix_init

  subroutine regmix_init(y,x,k,result,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k;type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call regmix_em(y,x,k,result,one_step_control(control),.true.)
  end subroutine regmix_init

  subroutine regmix_lambda_init(y,x,lambda_x,result,control)
    real(dp),intent(in)::y(:),x(:,:),lambda_x(:,:);type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call regmix_em_lambda(y,x,lambda_x,result,one_step_control(control),.true.)
  end subroutine regmix_lambda_init

  subroutine repnormmix_init(x,k,result,control)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k;type(mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call repnormmix_em(x,k,result,one_step_control(control))
  end subroutine repnormmix_init

  subroutine logisregmix_init(y,x,k,result,ntrials,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k;type(regression_mixture_result),intent(out)::result
    real(dp),intent(in),optional::ntrials(:);type(em_control),intent(in),optional::control
    call logisregmix_em(y,x,k,result,ntrials,one_step_control(control),.true.)
  end subroutine logisregmix_init

  subroutine poisregmix_init(y,x,k,result,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k;type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call poisregmix_em(y,x,k,result,one_step_control(control),.true.)
  end subroutine poisregmix_init

  subroutine segregmix_init(y,x,segmented_column,k,breakpoints,result,control)
    real(dp),intent(in)::y(:),x(:,:),breakpoints(:);integer,intent(in)::segmented_column,k
    type(regression_mixture_result),intent(out)::result;type(em_control),intent(in),optional::control
    call segregmix_em(y,x,segmented_column,k,breakpoints,result,one_step_control(control),.true.)
  end subroutine segregmix_init

  subroutine flaremix_init(y,x,k,result,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k;type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    call flaremix_em(y,x,k,result,one_step_control(control))
  end subroutine flaremix_init

  subroutine regmix_mixed_init(y,x,groups,k,result,random_effects,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::groups(:),k
    type(regression_mixture_result),intent(out)::result;real(dp),allocatable,intent(out)::random_effects(:,:)
    type(em_control),intent(in),optional::control
    call regmix_em_mixed(y,x,groups,k,result,random_effects,one_step_control(control),.true.)
  end subroutine regmix_mixed_init

  subroutine try_flare(y,x,k,result,control,nu)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k;type(regression_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control;real(dp),intent(in),optional::nu
    call flaremix_em(y,x,k,result,control,nu)
  end subroutine try_flare

  subroutine post_beta(draws,probability,mean,lower,upper,status)
    real(dp),intent(in)::draws(:,:,:),probability
    real(dp),allocatable,intent(out)::mean(:,:),lower(:,:),upper(:,:)
    integer,intent(out)::status
    call posterior_beta_intervals(draws,probability,mean,lower,upper,status)
  end subroutine post_beta

  subroutine regcr(draws,probability,mean,lower,upper,status)
    real(dp),intent(in)::draws(:,:,:),probability
    real(dp),allocatable,intent(out)::mean(:,:),lower(:,:),upper(:,:)
    integer,intent(out)::status
    call posterior_beta_intervals(draws,probability,mean,lower,upper,status)
  end subroutine regcr

  subroutine test_equality_vector(x,k,which_test,statistic,p_value,control)
    real(dp),intent(in)::x(:)
    integer,intent(in)::k,which_test
    real(dp),intent(out)::statistic,p_value
    type(em_control),intent(in),optional::control
    call test_equality_normal(x,k,which_test,statistic,p_value,control)
  end subroutine test_equality_vector

  subroutine test_equality_matrix(y,x,k,which_test,statistic,p_value,control)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k,which_test
    real(dp),intent(out)::statistic,p_value
    type(em_control),intent(in),optional::control
    call test_equality_regression(y,x,k,which_test,statistic,p_value,control)
  end subroutine test_equality_matrix

end module mixtools_compat
