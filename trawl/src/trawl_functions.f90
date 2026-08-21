module trawl_functions
  use trawl_kinds, only : dp
  use trawl_types, only : trawl_spec
  implicit none
  private
  public :: trawl_exp, trawl_dexp, trawl_supig, trawl_lm
  public :: acf_exp, acf_dexp, acf_supig, acf_lm
  public :: eval_trawl, trawl_tail_area
contains
  elemental real(dp) function trawl_exp(x,lambda) result(r)
    real(dp), intent(in) :: x,lambda
    r=exp(x*lambda)
  end function
  elemental real(dp) function trawl_dexp(x,w,lambda1,lambda2) result(r)
    real(dp), intent(in) :: x,w,lambda1,lambda2
    r=w*exp(lambda1*x)+(1.0_dp-w)*exp(lambda2*x)
  end function
  elemental real(dp) function trawl_supig(x,delta,gamma) result(r)
    real(dp), intent(in) :: x,delta,gamma
    real(dp) :: q
    q=1.0_dp-2.0_dp*x/gamma**2
    r=q**(-0.5_dp)*exp(delta*gamma*(1.0_dp-sqrt(q)))
  end function
  elemental real(dp) function trawl_lm(x,alpha,h) result(r)
    real(dp), intent(in) :: x,alpha,h
    r=(1.0_dp-x/alpha)**(-h)
  end function
  elemental real(dp) function acf_exp(x,lambda) result(r)
    real(dp), intent(in) :: x,lambda
    r=exp(-lambda*x)
  end function
  elemental real(dp) function acf_dexp(x,w,lambda1,lambda2) result(r)
    real(dp), intent(in) :: x,w,lambda1,lambda2
    r=(w*exp(-lambda1*x)/lambda1+(1.0_dp-w)*exp(-lambda2*x)/lambda2)/ &
      (w/lambda1+(1.0_dp-w)/lambda2)
  end function
  elemental real(dp) function acf_supig(x,delta,gamma) result(r)
    real(dp), intent(in) :: x,delta,gamma
    r=exp(delta*gamma*(1.0_dp-sqrt(1.0_dp+2.0_dp*x/gamma**2)))
  end function
  elemental real(dp) function acf_lm(x,alpha,h) result(r)
    real(dp), intent(in) :: x,alpha,h
    r=(1.0_dp+x/alpha)**(1.0_dp-h)
  end function
  real(dp) function eval_trawl(x,spec) result(r)
    real(dp), intent(in) :: x
    type(trawl_spec), intent(in) :: spec
    select case(trim(spec%kind))
    case('Exp','exp','EXP')
      r=trawl_exp(x,spec%lambda1)
    case('DExp','dexp','DEXP')
      r=trawl_dexp(x,spec%w,spec%lambda1,spec%lambda2)
    case('supIG','SUPIG','supig')
      r=trawl_supig(x,spec%delta,spec%gamma)
    case('LM','lm')
      r=trawl_lm(x,spec%alpha,spec%h)
    case default
      r=0.0_dp
    end select
  end function
  real(dp) function trawl_tail_area(z,spec) result(a)
    real(dp), intent(in) :: z
    type(trawl_spec), intent(in) :: spec
    real(dp) :: q
    select case(trim(spec%kind))
    case('Exp','exp','EXP')
      a=exp(spec%lambda1*z)/spec%lambda1
    case('DExp','dexp','DEXP')
      a=spec%w*exp(spec%lambda1*z)/spec%lambda1 + &
        (1.0_dp-spec%w)*exp(spec%lambda2*z)/spec%lambda2
    case('supIG','SUPIG','supig')
      if(spec%delta<=0.0_dp) then
        a=huge(1.0_dp)
      else
        q=sqrt(1.0_dp-2.0_dp*z/spec%gamma**2)
        a=spec%gamma/spec%delta*exp(spec%delta*spec%gamma*(1.0_dp-q))
      end if
    case('LM','lm')
      a=spec%alpha/(spec%h-1.0_dp)*(1.0_dp-z/spec%alpha)**(1.0_dp-spec%h)
    case default
      a=0.0_dp
    end select
  end function
end module trawl_functions
