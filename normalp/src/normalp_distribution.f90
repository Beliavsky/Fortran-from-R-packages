module normalp_distribution
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use normalp_special, only: dp, reg_gamma_p, gamma_quantile
  implicit none
  private
  public :: dnormp, pnormp, qnormp, rnormp
  interface dnormp
    module procedure dnormp_scalar, dnormp_vec
  end interface
  interface pnormp
    module procedure pnormp_scalar, pnormp_vec
  end interface
  interface qnormp
    module procedure qnormp_scalar, qnormp_vec
  end interface
contains
  pure function dnormp_scalar(x, mu, sigmap, p, log_pdf) result(y)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigmap, p
    logical, intent(in), optional :: log_pdf
    real(dp) :: y, muv, sv, pv, logy
    muv = 0.0_dp; if (present(mu)) muv = mu
    sv = 1.0_dp; if (present(sigmap)) sv = sigmap
    pv = 2.0_dp; if (present(p)) pv = p
    if (pv < 1.0_dp .or. sv <= 0.0_dp) then
      y = ieee_value(y, ieee_quiet_nan)
      return
    end if
    logy = -log(2.0_dp) - log(pv)/pv - log_gamma(1.0_dp + 1.0_dp/pv) - log(sv) &
           - abs(x-muv)**pv/(pv*sv**pv)
    if (present(log_pdf)) then
      if (log_pdf) then
        y = logy
      else
        y = exp(logy)
      end if
    else
      y = exp(logy)
    end if
  end function dnormp_scalar

  pure function dnormp_vec(x, mu, sigmap, p, log_pdf) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: mu, sigmap, p
    logical, intent(in), optional :: log_pdf
    real(dp) :: y(size(x)); integer :: i
    do i=1,size(x); y(i)=dnormp_scalar(x(i),mu,sigmap,p,log_pdf); end do
  end function dnormp_vec

  pure function pnormp_scalar(q, mu, sigmap, p, lower_tail, log_p) result(y)
    real(dp), intent(in) :: q
    real(dp), intent(in), optional :: mu, sigmap, p
    logical, intent(in), optional :: lower_tail, log_p
    real(dp) :: y, muv, sv, pv, z, zp
    logical :: lower, lp
    muv=0.0_dp; if(present(mu)) muv=mu
    sv=1.0_dp; if(present(sigmap)) sv=sigmap
    pv=2.0_dp; if(present(p)) pv=p
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    if (pv<1.0_dp .or. sv<=0.0_dp) then; y=ieee_value(y, ieee_quiet_nan); return; end if
    z=(q-muv)/sv
    zp=0.5_dp*reg_gamma_p(1.0_dp/pv,abs(z)**pv/pv)
    if(z<0.0_dp) then; y=0.5_dp-zp; else; y=0.5_dp+zp; end if
    if(.not.lower) y=1.0_dp-y
    if(lp) y=log(y)
  end function pnormp_scalar

  pure function pnormp_vec(q, mu, sigmap, p, lower_tail, log_p) result(y)
    real(dp),intent(in)::q(:); real(dp),intent(in),optional::mu,sigmap,p
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::y(size(q)); integer::i
    do i=1,size(q); y(i)=pnormp_scalar(q(i),mu,sigmap,p,lower_tail,log_p); end do
  end function pnormp_vec

  function qnormp_scalar(prob, mu, sigmap, p, lower_tail, log_p) result(q)
    real(dp), intent(in) :: prob
    real(dp), intent(in), optional :: mu, sigmap, p
    logical, intent(in), optional :: lower_tail, log_p
    real(dp) :: q, pr, muv, sv, pv, zp, qg, z
    logical :: lower, lp
    muv=0.0_dp; if(present(mu)) muv=mu
    sv=1.0_dp; if(present(sigmap)) sv=sigmap
    pv=2.0_dp; if(present(p)) pv=p
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    pr=prob
    if(lp) pr=exp(pr)
    if(.not.lower) pr=1.0_dp-pr
    if(pr<=0.0_dp) then; q=-huge(1.0_dp); return; end if
    if(pr>=1.0_dp) then; q=huge(1.0_dp); return; end if
    zp=2.0_dp*abs(pr-0.5_dp)
    if(zp==0.0_dp) then; q=muv; return; end if
    qg=gamma_quantile(zp,1.0_dp/pv,pv)
    z=qg**(1.0_dp/pv)
    if(pr<0.5_dp) z=-z
    q=muv+sv*z
  end function qnormp_scalar

  function qnormp_vec(prob, mu, sigmap, p, lower_tail, log_p) result(q)
    real(dp),intent(in)::prob(:); real(dp),intent(in),optional::mu,sigmap,p
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::q(size(prob)); integer::i
    do i=1,size(prob); q(i)=qnormp_scalar(prob(i),mu,sigmap,p,lower_tail,log_p); end do
  end function qnormp_vec

  subroutine rnormp(x, mu, sigmap, p)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: mu, sigmap, p
    real(dp) :: muv, sv, pv, u
    integer :: i
    muv=0.0_dp; if(present(mu)) muv=mu
    sv=1.0_dp; if(present(sigmap)) sv=sigmap
    pv=2.0_dp; if(present(p)) pv=p
    do i=1,size(x)
      call random_number(u)
      x(i)=qnormp_scalar(u,muv,sv,pv)
    end do
  end subroutine rnormp
end module normalp_distribution
