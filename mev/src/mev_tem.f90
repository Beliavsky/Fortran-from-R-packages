module mev_tem
  use mev_kinds, only: dp
  use mev_univariate, only: mev_fit_result, gpd_fit, gev_fit, gpd_ll, gev_ll, &
    gpd_infomat, gev_infomat
  use mev_profile, only: profile_result, gpd_profile, gev_profile
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: tem_profile_result, gpd_tem_profile, gev_tem_profile
  public :: gpd_vfun, gev_vfun, gpd_phi, gev_phi, gpd_dphi, gev_dphi

  type :: tem_profile_result
    real(dp), allocatable :: psi(:)
    real(dp), allocatable :: profile_ll(:)
    real(dp), allocatable :: r(:)
    real(dp), allocatable :: q(:)
    real(dp), allocatable :: rstar(:)
    real(dp), allocatable :: modified_tem_ll(:)
    real(dp) :: max_ll = -huge(1.0_dp)
    real(dp) :: mle_interest = 0.0_dp
    integer :: interest_index = 0
    integer :: convergence = 1
  end type tem_profile_result
contains

  subroutine gpd_vfun(par,dat,v)
    real(dp),intent(in)::par(2),dat(:)
    real(dp),allocatable,intent(out)::v(:,:)
    real(dp)::z,t
    integer::i
    allocate(v(size(dat),2));v=0.0_dp
    do i=1,size(dat)
      z=dat(i)/par(1);v(i,1)=z
      if(abs(par(2))<1.0e-6_dp) then
        v(i,2)=0.5_dp*par(1)*z*z
      else
        t=1.0_dp+par(2)*z
        if(t<=0.0_dp) then;v(i,2)=ieee_value(0.0_dp,ieee_quiet_nan)
        else;v(i,2)=par(1)*(t*log(t)-(t-1.0_dp))/(par(2)*par(2));end if
      end if
    end do
  end subroutine gpd_vfun

  subroutine gev_vfun(par,dat,v)
    real(dp),intent(in)::par(3),dat(:)
    real(dp),allocatable,intent(out)::v(:,:)
    real(dp)::z,t
    integer::i
    allocate(v(size(dat),3));v=0.0_dp
    do i=1,size(dat)
      z=(dat(i)-par(1))/par(2);v(i,1)=1.0_dp;v(i,2)=z
      if(abs(par(3))<1.0e-6_dp) then
        v(i,3)=0.5_dp*par(2)*z*z
      else
        t=1.0_dp+par(3)*z
        if(t<=0.0_dp) then;v(i,3)=ieee_value(0.0_dp,ieee_quiet_nan)
        else;v(i,3)=par(2)*(t*log(t)-(t-1.0_dp))/(par(3)*par(3));end if
      end if
    end do
  end subroutine gev_vfun

  subroutine gpd_phi(par,dat,v,phi)
    real(dp),intent(in)::par(2),dat(:),v(:,:)
    real(dp),intent(out)::phi(2)
    real(dp)::dx
    integer::i
    phi=0.0_dp
    if(par(1)<=0.0_dp.or.size(v,1)/=size(dat).or.size(v,2)/=2) then
      phi=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    do i=1,size(dat)
      if(par(1)+par(2)*dat(i)<=0.0_dp) then;phi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      dx=-(1.0_dp+par(2))/(par(1)+par(2)*dat(i))
      phi=phi+dx*v(i,:)
    end do
  end subroutine gpd_phi

  subroutine gev_phi(par,dat,v,phi)
    real(dp),intent(in)::par(3),dat(:),v(:,:)
    real(dp),intent(out)::phi(3)
    real(dp)::z,t,a,dx
    integer::i
    phi=0.0_dp
    if(par(2)<=0.0_dp.or.size(v,1)/=size(dat).or.size(v,2)/=3) then
      phi=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    do i=1,size(dat)
      z=(dat(i)-par(1))/par(2)
      if(abs(par(3))<1.0e-7_dp) then
        dx=(exp(-z)-1.0_dp)/par(2)
      else
        t=1.0_dp+par(3)*z
        if(t<=0.0_dp) then;phi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
        a=t**(-1.0_dp/par(3));dx=(a-1.0_dp-par(3))/(par(2)*t)
      end if
      phi=phi+dx*v(i,:)
    end do
  end subroutine gev_phi

  subroutine gpd_dphi(par,dat,v,dphi)
    real(dp),intent(in)::par(2),dat(:),v(:,:)
    real(dp),intent(out)::dphi(2,2)
    real(dp)::pp(2),pm(2),fp(2),fm(2),h
    integer::j,it
    logical::okp,okm
    do j=1,2
      if(j==1)then;h=5.0e-5_dp*max(par(1),1.0_dp);else;h=1.0e-4_dp*max(abs(par(2)),1.0_dp);end if
      do it=1,10
        pp=par;pm=par;pp(j)=pp(j)+h;pm(j)=pm(j)-h
        call gpd_phi(pp,dat,v,fp);call gpd_phi(pm,dat,v,fm)
        okp=all(ieee_is_finite(fp));okm=all(ieee_is_finite(fm))
        if(okp.and.okm)exit
        h=0.5_dp*h
      end do
      if(.not.(okp.and.okm))then;dphi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      dphi(:,j)=(fp-fm)/(2.0_dp*h)
    end do
  end subroutine gpd_dphi

  subroutine gev_dphi(par,dat,v,dphi)
    real(dp),intent(in)::par(3),dat(:),v(:,:)
    real(dp),intent(out)::dphi(3,3)
    real(dp)::pp(3),pm(3),fp(3),fm(3),h
    integer::j,it
    logical::okp,okm
    do j=1,3
      select case(j)
      case(1);h=5.0e-5_dp*max(par(2),max(abs(par(1)),1.0_dp))
      case(2);h=5.0e-5_dp*max(par(2),1.0_dp)
      case default;h=1.0e-4_dp*max(abs(par(3)),1.0_dp)
      end select
      do it=1,10
        pp=par;pm=par;pp(j)=pp(j)+h;pm(j)=pm(j)-h
        call gev_phi(pp,dat,v,fp);call gev_phi(pm,dat,v,fm)
        okp=all(ieee_is_finite(fp));okm=all(ieee_is_finite(fm))
        if(okp.and.okm)exit
        h=0.5_dp*h
      end do
      if(.not.(okp.and.okm))then;dphi=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      dphi(:,j)=(fp-fm)/(2.0_dp*h)
    end do
  end subroutine gev_dphi

  subroutine gpd_tem_profile(dat,psi,param,res)
    real(dp),intent(in)::dat(:),psi(:)
    character(len=*),intent(in)::param
    type(tem_profile_result),intent(out)::res
    type(mev_fit_result)::fit
    type(profile_result)::pr
    real(dp),allocatable::v(:,:)
    real(dp)::mle(2),par(2),phimle(2),ph(2),dphimle(2,2),dph(2,2),info(2,2)
    real(dp)::infomle(2,2),qnum,logq,qmle,rval,detn,detfull
    integer::ind,i,other
    call gpd_fit(dat,fit);mle=fit%estimate;call gpd_profile(dat,psi,param,pr)
    if(trim(adjustl(param))=='scale')then;ind=1;else if(trim(adjustl(param))=='shape')then;ind=2;else;return;end if
    allocate(res%psi(size(psi)),res%profile_ll(size(psi)),res%r(size(psi)),res%q(size(psi)), &
      res%rstar(size(psi)),res%modified_tem_ll(size(psi)))
    res%psi=psi;res%profile_ll=pr%loglik;res%max_ll=fit%loglik;res%mle_interest=mle(ind);res%interest_index=ind
    call gpd_vfun(mle,dat,v);call gpd_phi(mle,dat,v,phimle);call gpd_dphi(mle,dat,v,dphimle)
    call gpd_infomat(mle,dat,infomle,.false.)
    detfull=det2(dphimle);qmle=-log(max(abs(detfull),tiny(1.0_dp)))+0.5_dp*log(max(det2(infomle),tiny(1.0_dp)))
    other=3-ind;res%convergence=0
    do i=1,size(psi)
      if(ind==1)then;par=[psi(i),pr%nuisance(i,1)];else;par=[pr%nuisance(i,1),psi(i)];end if
      call gpd_phi(par,dat,v,ph);call gpd_dphi(par,dat,v,dph);call gpd_infomat(par,dat,info,.false.)
      qnum=phimle(ind)-ph(ind)
      ! determinant after replacing the interest row by phi(mle)-phi(par)
      if(ind==1) qnum=(phimle(1)-ph(1))*dph(2,2)-(phimle(2)-ph(2))*dph(2,1)
      if(ind==2) qnum=-((phimle(1)-ph(1))*dph(1,2)-(phimle(2)-ph(2))*dph(1,1))
      detn=max(info(other,other),tiny(1.0_dp))
      logq=-0.5_dp*log(detn)+log(max(abs(qnum),tiny(1.0_dp)))+qmle
      rval=sign(1.0_dp,mle(ind)-psi(i))*sqrt(max(0.0_dp,2.0_dp*(fit%loglik-pr%loglik(i))))
      res%r(i)=rval;res%q(i)=sign(exp(logq),qnum)
      if(abs(rval)>1.0e-10_dp)then;res%rstar(i)=rval+(logq-log(abs(rval)))/rval;else;res%rstar(i)=0.0_dp;end if
      res%modified_tem_ll(i)=pr%loglik(i)+0.5_dp*log(detn)-log(max(abs(dph(other,other)),tiny(1.0_dp)))
    end do
  end subroutine gpd_tem_profile

  subroutine gev_tem_profile(dat,psi,param,res)
    real(dp),intent(in)::dat(:),psi(:)
    character(len=*),intent(in)::param
    type(tem_profile_result),intent(out)::res
    type(mev_fit_result)::fit
    type(profile_result)::pr
    real(dp),allocatable::v(:,:)
    real(dp)::mle(3),par(3),phimle(3),ph(3),dphimle(3,3),dph(3,3),info(3,3),infomle(3,3)
    real(dp)::row(3),qnum,logq,qmle,rval,detn,detfull,detdn
    integer::ind,i,nidx(2)
    character(len=16)::pa
    pa=trim(adjustl(param));select case(pa);case('loc');ind=1;case('scale');ind=2;case('shape');ind=3;case default;return;end select
    call gev_fit(dat,fit);mle=fit%estimate;call gev_profile(dat,psi,param,pr)
    allocate(res%psi(size(psi)),res%profile_ll(size(psi)),res%r(size(psi)),res%q(size(psi)), &
      res%rstar(size(psi)),res%modified_tem_ll(size(psi)))
    res%psi=psi;res%profile_ll=pr%loglik;res%max_ll=fit%loglik;res%mle_interest=mle(ind);res%interest_index=ind
    call gev_vfun(mle,dat,v);call gev_phi(mle,dat,v,phimle);call gev_dphi(mle,dat,v,dphimle)
    call gev_infomat(mle,dat,infomle,.false.);detfull=det3(dphimle)
    qmle=-log(max(abs(detfull),tiny(1.0_dp)))+0.5_dp*log(max(det3(infomle),tiny(1.0_dp)))
    call nuisance_indices(ind,nidx);res%convergence=0
    do i=1,size(psi)
      select case(ind)
      case(1);par=[psi(i),pr%nuisance(i,1),pr%nuisance(i,2)]
      case(2);par=[pr%nuisance(i,1),psi(i),pr%nuisance(i,2)]
      case(3);par=[pr%nuisance(i,1),pr%nuisance(i,2),psi(i)]
      end select
      call gev_phi(par,dat,v,ph);call gev_dphi(par,dat,v,dph);call gev_infomat(par,dat,info,.false.)
      row=phimle-ph
      qnum=oriented_det3(row,dph,nidx,ind)
      detn=det2_sub(info,nidx);detdn=det2_rows_cols(dph,nidx,nidx)
      logq=-0.5_dp*log(max(detn,tiny(1.0_dp)))+log(max(abs(qnum),tiny(1.0_dp)))+qmle
      rval=sign(1.0_dp,mle(ind)-psi(i))*sqrt(max(0.0_dp,2.0_dp*(fit%loglik-pr%loglik(i))))
      res%r(i)=rval;res%q(i)=sign(exp(logq),qnum)
      if(abs(rval)>1.0e-10_dp)then;res%rstar(i)=rval+(logq-log(abs(rval)))/rval;else;res%rstar(i)=0.0_dp;end if
      res%modified_tem_ll(i)=pr%loglik(i)+0.5_dp*log(max(detn,tiny(1.0_dp)))-log(max(abs(detdn),tiny(1.0_dp)))
    end do
  end subroutine gev_tem_profile

  pure subroutine nuisance_indices(ind,idx)
    integer,intent(in)::ind
    integer,intent(out)::idx(2)
    select case(ind);case(1);idx=[2,3];case(2);idx=[1,3];case default;idx=[1,2];end select
  end subroutine
  pure real(dp) function det2(a) result(d);real(dp),intent(in)::a(2,2);d=a(1,1)*a(2,2)-a(1,2)*a(2,1);end function
  pure real(dp) function det3(a) result(d)
    real(dp),intent(in)::a(3,3)
    d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+ &
      a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function
  pure real(dp) function det2_sub(a,idx) result(d)
    real(dp),intent(in)::a(3,3);integer,intent(in)::idx(2)
    d=a(idx(1),idx(1))*a(idx(2),idx(2))-a(idx(1),idx(2))*a(idx(2),idx(1))
  end function
  pure real(dp) function det2_rows_cols(a,rows,cols) result(d)
    real(dp),intent(in)::a(3,3);integer,intent(in)::rows(2),cols(2)
    d=a(rows(1),cols(1))*a(rows(2),cols(2))-a(rows(1),cols(2))*a(rows(2),cols(1))
  end function
  pure real(dp) function oriented_det3(row,dphi,nidx,ind) result(d)
    real(dp),intent(in)::row(3),dphi(3,3);integer,intent(in)::nidx(2),ind
    real(dp)::a(3,3)
    a(1,:)=row;a(2,:)=dphi(nidx(1),:);a(3,:)=dphi(nidx(2),:)
    d=det3(a);if(mod(ind,2)==0)d=-d
  end function
end module mev_tem
