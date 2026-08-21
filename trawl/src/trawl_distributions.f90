module trawl_distributions
  use trawl_kinds, only : dp
  use trawl_rng, only : rnbinom_scalar,rlogarithmic_scalar,rbinom1_scalar
  use trawl_types, only : trawl_ok,trawl_invalid_argument
  implicit none
  private
  public :: bivariate_nbsim,bivariate_lsdsim,trivariate_lsdsim
  public :: lsd_mean,lsd_var,modlsd_mean,modlsd_var
  public :: bivlsd_cor,bivlsd_cov,bivmodlsd_cov,bivmodlsd_cor
contains
  pure real(dp) function lsd_mean(p) result(r)
    real(dp),intent(in)::p
    r=p/((1.0_dp-p)*(-log(1.0_dp-p)))
  end function
  pure real(dp) function lsd_var(p) result(r)
    real(dp),intent(in)::p
    r=(-p)*(p+log(1.0_dp-p))/((1.0_dp-p)*log(1.0_dp-p))**2
  end function
  pure real(dp) function modlsd_mean(delta,p) result(r)
    real(dp),intent(in)::delta,p
    r=(1.0_dp-delta)*lsd_mean(p)
  end function
  pure real(dp) function modlsd_var(delta,p) result(r)
    real(dp),intent(in)::delta,p
    real(dp)::eb,eb2,varb
    eb=1.0_dp-delta
    eb2=delta*(1.0_dp-delta)+(1.0_dp-delta)**2
    varb=delta*(1.0_dp-delta)
    r=eb2*lsd_var(p)+varb*lsd_mean(p)**2
  end function
  pure real(dp) function bivlsd_cor(p1,p2) result(r)
    real(dp),intent(in)::p1,p2
    real(dp)::gamma,l
    gamma=1.0_dp-p1-p2; l=-log(gamma)
    r=sign(1.0_dp,l-1.0_dp)*sqrt(p1*p2*(l-1.0_dp)**2/ &
      ((p1*(l-1.0_dp)+gamma*l)*(p2*(l-1.0_dp)+gamma*l)))
  end function
  pure real(dp) function bivlsd_cov(p1,p2) result(r)
    real(dp),intent(in)::p1,p2
    real(dp)::tp1,tp2,d1,d2
    tp1=p1/(1.0_dp-p2); tp2=p2/(1.0_dp-p1)
    d1=log(1.0_dp-p2)/log(1.0_dp-p1-p2)
    d2=log(1.0_dp-p1)/log(1.0_dp-p1-p2)
    r=bivlsd_cor(p1,p2)*sqrt(modlsd_var(d1,tp1)*modlsd_var(d2,tp2))
  end function
  pure real(dp) function bivmodlsd_cov(delta,p1,p2) result(r)
    real(dp),intent(in)::delta,p1,p2
    real(dp)::eb2,varb,tp1,tp2,d1,d2,ea1,ea2
    eb2=delta*(1.0_dp-delta)+(1.0_dp-delta)**2
    varb=delta*(1.0_dp-delta)
    tp1=p1/(1.0_dp-p2); tp2=p2/(1.0_dp-p1)
    d1=log(1.0_dp-p2)/log(1.0_dp-p1-p2)
    d2=log(1.0_dp-p1)/log(1.0_dp-p1-p2)
    ea1=(1.0_dp-d1)*lsd_mean(tp1); ea2=(1.0_dp-d2)*lsd_mean(tp2)
    r=eb2*bivlsd_cov(p1,p2)+ea1*ea2*varb
  end function
  pure real(dp) function bivmodlsd_cor(delta,p1,p2) result(r)
    real(dp),intent(in)::delta,p1,p2
    real(dp)::tp1,tp2,d1,d2
    tp1=p1/(1.0_dp-p2); tp2=p2/(1.0_dp-p1)
    d1=log(1.0_dp-p2)/log(1.0_dp-p1-p2)
    d2=log(1.0_dp-p1)/log(1.0_dp-p1-p2)
    r=bivmodlsd_cov(delta,p1,p2)/sqrt(modlsd_var(d1,tp1)*modlsd_var(d2,tp2))
  end function

  subroutine bivariate_nbsim(n,kappa,p1,p2,x,status)
    integer,intent(in)::n
    real(dp),intent(in)::kappa,p1,p2
    integer,allocatable,intent(out)::x(:,:)
    integer,intent(out),optional::status
    integer::i,s
    s=trawl_ok
    if(n<=0 .or. kappa<=0.0_dp .or. p1<=0.0_dp .or. p2<=0.0_dp .or. &
       p1>=1.0_dp .or. p2>=1.0_dp .or. p1+p2>=1.0_dp) s=trawl_invalid_argument
    if(s/=trawl_ok) then
      allocate(x(0,2)); if(present(status)) status=s; return
    end if
    allocate(x(n,2))
    do i=1,n
      x(i,1)=rnbinom_scalar(kappa,1.0_dp-p1/(1.0_dp-p2))
      x(i,2)=rnbinom_scalar(kappa+real(x(i,1),dp),1.0_dp-p2)
    end do
    if(present(status)) status=s
  end subroutine

  subroutine bivariate_lsdsim(n,p1,p2,x,status)
    integer,intent(in)::n
    real(dp),intent(in)::p1,p2
    integer,allocatable,intent(out)::x(:,:)
    integer,intent(out),optional::status
    integer::i,s,b
    real(dp)::delta1,tp1
    s=trawl_ok
    if(n<=0 .or. p1<=0.0_dp .or. p2<=0.0_dp .or. p1>=1.0_dp .or. p2>=1.0_dp .or. p1+p2>=1.0_dp) s=trawl_invalid_argument
    if(s/=trawl_ok) then
      allocate(x(0,2)); if(present(status)) status=s; return
    end if
    delta1=log(1.0_dp-p2)/log(1.0_dp-p1-p2); tp1=p1/(1.0_dp-p2)
    allocate(x(n,2))
    do i=1,n
      b=rbinom1_scalar(1.0_dp-delta1)
      x(i,1)=rlogarithmic_scalar(tp1)*b
      if(x(i,1)==0) then
        x(i,2)=rlogarithmic_scalar(p2)
      else
        x(i,2)=rnbinom_scalar(real(x(i,1),dp),1.0_dp-p2)
      end if
    end do
    if(present(status)) status=s
  end subroutine

  subroutine trivariate_lsdsim(n,p1,p2,p3,x,status)
    integer,intent(in)::n
    real(dp),intent(in)::p1,p2,p3
    integer,allocatable,intent(out)::x(:,:)
    integer,intent(out),optional::status
    integer::i,s,b
    real(dp)::delta1,tp1
    integer,allocatable::y(:,:)
    s=trawl_ok
    if(n<=0 .or. min(p1,p2,p3)<=0.0_dp .or. max(p1,p2,p3)>=1.0_dp .or. p1+p2+p3>=1.0_dp) s=trawl_invalid_argument
    if(s/=trawl_ok) then
      allocate(x(0,3)); if(present(status)) status=s; return
    end if
    delta1=log(1.0_dp-p2-p3)/log(1.0_dp-p1-p2-p3); tp1=p1/(1.0_dp-p2-p3)
    allocate(x(n,3))
    do i=1,n
      b=rbinom1_scalar(1.0_dp-delta1)
      x(i,1)=rlogarithmic_scalar(tp1)*b
      if(x(i,1)==0) then
        call bivariate_lsdsim(1,p2,p3,y,s)
      else
        call bivariate_nbsim(1,real(x(i,1),dp),p2,p3,y,s)
      end if
      x(i,2:3)=y(1,:)
    end do
    if(present(status)) status=s
  end subroutine
end module trawl_distributions
