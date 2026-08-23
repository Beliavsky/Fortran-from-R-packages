module mev_reparam
  use mev_kinds, only: dp
  use mev_math, only: digamma_mev, finite_diff_gradient, finite_diff_hessian
  use mev_univariate, only: gpd_ll, gev_ll
  implicit none
  private
  public :: gpde_ll, gpde_score, gpde_infomat, gpdr_ll, gpdr_score, gpdr_infomat
  public :: gevr_ll, gevr_score, gevr_infomat, gpdn_ll, gpdn_score, gpdn_infomat
  public :: gevn_ll, gevn_score, gevn_infomat
  public :: rrlarg, rlarg_ll, rlarg_score, rlarg_infomat
  public :: pp_ll, pp_score, pp_infomat
contains
  pure real(dp) function gpde_scale(es,xi,m) result(sig)
    real(dp),intent(in)::es,xi,m
    if(abs(xi)>1.0e-9_dp)then
      sig=es*(1.0_dp-xi)*xi/(m**xi-1.0_dp+xi)
    else
      sig=es/(log(m)+1.0_dp)
    end if
  end function gpde_scale

  real(dp) function gpde_ll(par,dat,m) result(ll)
    real(dp),intent(in)::par(:),dat(:),m
    real(dp)::sig
    if(size(par)/=2.or.par(1)<=0.0_dp.or.par(2)>1.0_dp.or.m<=1.0_dp)then;ll=-huge(1.0_dp);return;end if
    sig=gpde_scale(par(1),par(2),m)
    if(sig<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    ll=gpd_ll([sig,par(2)],dat)
  end function gpde_ll

  subroutine gpde_score(par,dat,m,score)
    real(dp),intent(in)::par(:),dat(:),m
    real(dp),intent(out)::score(2)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpde_ll(x,dat,m);end function
  end subroutine gpde_score

  subroutine gpde_infomat(par,dat,m,info)
    real(dp),intent(in)::par(:),dat(:),m
    real(dp),intent(out)::info(2,2)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpde_ll(x,dat,m);end function
  end subroutine gpde_infomat

  pure real(dp) function gpdr_scale(ym,xi,m) result(sig)
    real(dp),intent(in)::ym,xi,m
    if(abs(xi)>1.0e-8_dp)then;sig=ym*xi/(m**xi-1.0_dp);else;sig=ym/log(m);end if
  end function gpdr_scale

  real(dp) function gpdr_ll(par,dat,m) result(ll)
    real(dp),intent(in)::par(:),dat(:),m
    real(dp)::sig
    if(size(par)/=2.or.par(1)<=0.0_dp.or.m<=1.0_dp)then;ll=-huge(1.0_dp);return;end if
    sig=gpdr_scale(par(1),par(2),m);if(sig<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    ll=gpd_ll([sig,par(2)],dat)
  end function gpdr_ll

  subroutine gpdr_score(par,dat,m,score)
    real(dp),intent(in)::par(:),dat(:),m;real(dp),intent(out)::score(2)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpdr_ll(x,dat,m);end function
  end subroutine gpdr_score

  subroutine gpdr_infomat(par,dat,m,info)
    real(dp),intent(in)::par(:),dat(:),m;real(dp),intent(out)::info(2,2)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpdr_ll(x,dat,m);end function
  end subroutine gpdr_infomat

  real(dp) function gevr_ll(par,dat,p) result(ll)
    real(dp),intent(in)::par(:),dat(:),p
    real(dp)::mu,sig,xi,z
    if(size(par)/=3.or.par(2)<=0.0_dp.or.p<=0.0_dp.or.p>=1.0_dp)then;ll=-huge(1.0_dp);return;end if
    z=par(1);sig=par(2);xi=par(3)
    if(abs(xi)>1.0e-8_dp)then
      mu=z+sig/xi*(1.0_dp-(-log(1.0_dp-p))**(-xi))
    else
      mu=z+sig*log(-log(1.0_dp-p))
    end if
    ll=gev_ll([mu,sig,xi],dat)
  end function gevr_ll

  subroutine gevr_score(par,dat,p,score)
    real(dp),intent(in)::par(:),dat(:),p;real(dp),intent(out)::score(3)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gevr_ll(x,dat,p);end function
  end subroutine gevr_score

  subroutine gevr_infomat(par,dat,p,info)
    real(dp),intent(in)::par(:),dat(:),p;real(dp),intent(out)::info(3,3)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gevr_ll(x,dat,p);end function
  end subroutine gevr_infomat

  pure real(dp) function gpdn_scale(z,xi,nblock) result(sig)
    real(dp),intent(in)::z,xi,nblock
    real(dp),parameter::euler=0.57721566490153231044_dp
    if(abs(xi)>1.0e-8_dp)then
      sig=z*xi/(exp(log_gamma(nblock+1.0_dp)+log_gamma(1.0_dp-xi)-log_gamma(nblock-xi+1.0_dp))-1.0_dp)
    else
      sig=z/(euler+digamma_mev(nblock+1.0_dp))
    end if
  end function gpdn_scale

  real(dp) function gpdn_ll(par,dat,nblock) result(ll)
    real(dp),intent(in)::par(:),dat(:),nblock
    real(dp)::sig
    if(size(par)/=2.or.nblock<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    sig=gpdn_scale(par(1),par(2),nblock);if(sig<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    ll=gpd_ll([sig,par(2)],dat)
  end function gpdn_ll

  subroutine gpdn_score(par,dat,nblock,score)
    real(dp),intent(in)::par(:),dat(:),nblock;real(dp),intent(out)::score(2)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpdn_ll(x,dat,nblock);end function
  end subroutine gpdn_score

  subroutine gpdn_infomat(par,dat,nblock,info)
    real(dp),intent(in)::par(:),dat(:),nblock;real(dp),intent(out)::info(2,2)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gpdn_ll(x,dat,nblock);end function
  end subroutine gpdn_infomat

  real(dp) function gevn_ll(par,dat,nblock,q,qty) result(ll)
    real(dp),intent(in)::par(:),dat(:),nblock,q
    character(len=*),intent(in)::qty
    real(dp)::mu,z,xi,sig
    if(size(par)/=3.or.nblock<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    mu=par(1);z=par(2);xi=par(3)
    if(abs(xi)>1.0e-8_dp)then
      if(trim(qty)=='quantile')then
        sig=(z-mu)*xi/(nblock**xi*(log(1.0_dp/q))**(-xi)-1.0_dp)
      else
        sig=(z-mu)*xi/(nblock**xi*gamma(1.0_dp-xi)-1.0_dp)
      end if
    else
      if(trim(qty)=='quantile')then
        sig=(z-mu)/(log(nblock)-log(-log(q)))
      else
        sig=(z-mu)/(-digamma_mev(1.0_dp)+log(nblock))
      end if
    end if
    if(sig<=0.0_dp)then;ll=-huge(1.0_dp);else;ll=gev_ll([mu,sig,xi],dat);end if
  end function gevn_ll

  subroutine gevn_score(par,dat,nblock,q,qty,score)
    real(dp),intent(in)::par(:),dat(:),nblock,q;character(len=*),intent(in)::qty;real(dp),intent(out)::score(3)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gevn_ll(x,dat,nblock,q,qty);end function
  end subroutine gevn_score

  subroutine gevn_infomat(par,dat,nblock,q,qty,info)
    real(dp),intent(in)::par(:),dat(:),nblock,q;character(len=*),intent(in)::qty;real(dp),intent(out)::info(3,3)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=gevn_ll(x,dat,nblock,q,qty);end function
  end subroutine gevn_infomat

  subroutine rrlarg(n,r,loc,scale,shape,sample)
    integer,intent(in)::n,r;real(dp),intent(in)::loc,scale,shape
    real(dp),intent(out)::sample(n,r)
    real(dp)::cs,u
    integer::i,j
    do i=1,n
      cs=0.0_dp
      do j=1,r
        call random_number(u);u=max(u,tiny(1.0_dp));cs=cs-log(u)
        if(abs(shape)>1.0e-12_dp)then
          sample(i,j)=loc+scale*(cs**(-shape)-1.0_dp)/shape
        else
          sample(i,j)=loc-scale*log(cs)
        end if
      end do
    end do
  end subroutine rrlarg

  real(dp) function rlarg_ll(par,dat) result(ll)
    real(dp),intent(in)::par(:),dat(:,:)
    real(dp)::loc,scale,shape,t
    integer::i,n,r
    if(size(par)/=3.or.par(2)<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    loc=par(1);scale=par(2);shape=par(3);n=size(dat,1);r=size(dat,2);ll=0.0_dp
    do i=1,n
      if(any(dat(i,1:r-1)<dat(i,2:r)))then;ll=-huge(1.0_dp);return;end if
      if(abs(shape)>1.0e-7_dp)then
        t=1.0_dp+shape*(dat(i,r)-loc)/scale
        if(t<=0.0_dp.or.any(1.0_dp+shape*(dat(i,:)-loc)/scale<=0.0_dp))then;ll=-huge(1.0_dp);return;end if
        ll=ll-t**(-1.0_dp/shape)-real(r,dp)*log(scale) &
          -(1.0_dp/shape+1.0_dp)*sum(log(1.0_dp+shape*(dat(i,:)-loc)/scale))
      else
        ll=ll-exp((loc-dat(i,r))/scale)-real(r,dp)*log(scale)-sum((dat(i,:)-loc)/scale)
      end if
    end do
  end function rlarg_ll

  subroutine rlarg_score(par,dat,score)
    real(dp),intent(in)::par(:),dat(:,:);real(dp),intent(out)::score(3)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=rlarg_ll(x,dat);end function
  end subroutine rlarg_score

  subroutine rlarg_infomat(par,dat,info)
    real(dp),intent(in)::par(:),dat(:,:);real(dp),intent(out)::info(3,3)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=rlarg_ll(x,dat);end function
  end subroutine rlarg_infomat

  real(dp) function pp_ll(par,dat,u,np) result(ll)
    real(dp),intent(in)::par(:),dat(:),u,np
    real(dp)::loc,scale,shape,t
    if(size(par)/=3.or.par(2)<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
    loc=par(1);scale=par(2);shape=par(3)
    if(abs(shape)>1.0e-6_dp)then
      t=1.0_dp+shape*(u-loc)/scale
      if(t<=0.0_dp.or.any(1.0_dp+shape*(dat-loc)/scale<=0.0_dp))then;ll=-huge(1.0_dp);return;end if
      ll=-np*t**(-1.0_dp/shape)-real(size(dat),dp)*log(scale) &
        -(1.0_dp+1.0_dp/shape)*sum(log(1.0_dp+shape*(dat-loc)/scale))
    else
      ll=-np*exp((loc-u)/scale)-real(size(dat),dp)*log(scale)+sum(loc-dat)/scale
    end if
  end function pp_ll

  subroutine pp_score(par,dat,u,np,score)
    real(dp),intent(in)::par(:),dat(:),u,np;real(dp),intent(out)::score(3)
    call finite_diff_gradient(obj,par,score,1.0e-5_dp)
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=pp_ll(x,dat,u,np);end function
  end subroutine pp_score

  subroutine pp_infomat(par,dat,u,np,info)
    real(dp),intent(in)::par(:),dat(:),u,np;real(dp),intent(out)::info(3,3)
    call finite_diff_hessian(obj,par,info,1.0e-4_dp);info=-info
  contains
    function obj(x) result(v);real(dp),intent(in)::x(:);real(dp)::v;v=pp_ll(x,dat,u,np);end function
  end subroutine pp_infomat
end module mev_reparam
