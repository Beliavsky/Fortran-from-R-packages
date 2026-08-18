! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran
  use runuran_kinds
  use runuran_rng
  use runuran_math
  use runuran_distributions
  use runuran_methods
  use runuran_multivariate
  implicit none
  private
  public :: dp, i8, pi
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_exponential
  public :: rng_gamma, rng_beta, rng_poisson, rng_binomial, rng_geometric
  public :: rng_negative_binomial, rng_cauchy, rng_chisq, rng_student_t, rng_f
  public :: continuous_distribution, discrete_distribution, multivariate_distribution
  public :: unuran_generator, mv_generator
  public :: ud_continuous, ud_continuous_cdf, ud_continuous_logpdf, ud_discrete
  public :: ud_probability_vector, udnorm, udbeta, udcauchy, udchi, udchisq
  public :: udexp, udf, udfrechet, udgamma, udgumbel, udig, udlaplace, udlnorm
  public :: udlogis, udlomax, udpareto, udpowerexp, udrayleigh, udslash, udt
  public :: udtriang, udweibull, udburr, udgig, udgiga, udhyperbolic, udghyp
  public :: udvg, udmeixner, udplanck
  public :: udbinom, udgeom, udhyper, udlogarithmic, udnbinom, udpois, udzipf
  public :: udmultinormal, udmultistudent, udmulticauchy, udmultiexponential
  public :: udmultivariate
  public :: pinv_new, ars_new, arou_new, srou_new, tdr_new, itdr_new, tabl_new
  public :: dari_new, dau_new, dgt_new, mixt_new, unuran_new
  public :: hitro_new, vnrou_new, gibbs_new
  public :: METHOD_PINV, METHOD_ARS, METHOD_AROU, METHOD_SROU, METHOD_TDR
  public :: METHOD_ITDR, METHOD_TABL, METHOD_DARI, METHOD_DAU, METHOD_DGT
  public :: ud, up, uq, ur, ur_int, unuran_is_inversion

  interface ud
    module procedure ud_cont_value, ud_discr_value, ud_gen_value
  end interface
  interface up
    module procedure up_cont_value, up_discr_value, up_gen_value
  end interface
contains
  real(dp) function ud_cont_value(d,x,logscale) result(v)
    type(continuous_distribution),intent(in)::d
    real(dp),intent(in)::x
    logical,intent(in),optional::logscale
    logical::lg
    lg=.false.;if(present(logscale))lg=logscale
    if(lg)then;v=d%logpdf(x);else;v=d%pdf(x);end if
  end function
  real(dp) function ud_discr_value(d,k,logscale) result(v)
    type(discrete_distribution),intent(in)::d
    integer,intent(in)::k
    logical,intent(in),optional::logscale
    logical::lg
    lg=.false.;if(present(logscale))lg=logscale
    v=d%pmf(k);if(lg)v=log(max(v,tiny(1.0_dp)))
  end function
  real(dp) function up_cont_value(d,x) result(v)
    type(continuous_distribution),intent(in)::d;real(dp),intent(in)::x
    v=d%cdf(x)
  end function
  real(dp) function up_discr_value(d,k) result(v)
    type(discrete_distribution),intent(in)::d;integer,intent(in)::k
    v=d%cdf(k)
  end function
  real(dp) function ud_gen_value(g,x,logscale) result(v)
    type(unuran_generator),intent(in)::g
    real(dp),intent(in)::x
    logical,intent(in),optional::logscale
    logical::lg
    lg=.false.
    if(present(logscale))lg=logscale
    if(g%is_discrete)then
      v=g%discr%pmf(nint(x))
      if(lg)v=log(max(v,tiny(1.0_dp)))
    else
      if(lg)then
        v=g%cont%logpdf(x)
      else
        v=g%cont%pdf(x)
      end if
    end if
  end function ud_gen_value

  real(dp) function up_gen_value(g,x) result(v)
    type(unuran_generator),intent(in)::g
    real(dp),intent(in)::x
    if(g%is_discrete)then
      v=g%discr%cdf(floor(x))
    else
      v=g%cont%cdf(x)
    end if
  end function up_gen_value

  real(dp) function uq(g,p) result(x)
    type(unuran_generator),intent(in)::g;real(dp),intent(in)::p
    if(g%is_discrete)then;x=real(g%discr%quantile(p),dp);else;x=g%cont%quantile(p);end if
  end function
  real(dp) function ur(g,rng) result(x)
    type(unuran_generator),intent(inout)::g;type(rng_state),intent(inout)::rng
    x=g%sample(rng)
  end function
  integer function ur_int(g,rng) result(x)
    type(unuran_generator),intent(inout)::g;type(rng_state),intent(inout)::rng
    x=g%sample_int(rng)
  end function
  logical function unuran_is_inversion(g) result(tf)
    type(unuran_generator),intent(in)::g
    tf=(g%method==METHOD_PINV .or. g%method==METHOD_DARI .or. g%method==METHOD_DAU)
  end function unuran_is_inversion
end module runuran
