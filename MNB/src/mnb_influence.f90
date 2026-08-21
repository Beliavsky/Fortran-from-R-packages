! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_influence
  use mnb_kinds, only : dp
  use mnb_types, only : mnb_fit_result,mnb_global_result,mnb_local_result
  use mnb_core, only : fit_mnb,mnb_loglik,mnb_cluster_sums
  use mnb_math, only : invert_matrix,numerical_hessian,symmetric_eigen_jacobi,sample_sd_unique
  implicit none
  private
  integer,parameter,public :: local_weight=1,local_weight_obs=2,local_covariate=3,local_dispersion=4
  public :: global_influence_mnb, local_influence_mnb
contains
  function global_influence_mnb(start,y,x,n,mi,offset) result(out)
    real(dp),intent(in)::start(:),y(:),x(:,:)
    integer,intent(in)::n,mi
    real(dp),intent(in),optional::offset(:)
    type(mnb_global_result)::out
    type(mnb_fit_result)::base,del
    real(dp),allocatable::yd(:),xd(:,:),od(:),infoinv(:,:),dpar(:)
    logical::ok
    integer::i,l1,l2,nn,p
    nn=n*mi;p=size(x,2);base=fit_mnb(start,y,x,n,mi,offset)
    allocate(out%cook_distance(n),out%likelihood_displacement(n),infoinv(size(start),size(start)),dpar(size(start)))
    call invert_matrix(-base%hessian,infoinv,ok);if(.not.ok)error stop 'global_influence_mnb: singular information'
    do i=1,n
      allocate(yd(nn-mi),xd(nn-mi,p));l1=(i-1)*mi;l2=i*mi
      if(l1>0)then;yd(1:l1)=y(1:l1);xd(1:l1,:)=x(1:l1,:);end if
      if(l2<nn)then;yd(l1+1:)=y(l2+1:);xd(l1+1:,:)=x(l2+1:,:);end if
      if(present(offset))then
        allocate(od(nn-mi));if(l1>0)od(1:l1)=offset(1:l1);if(l2<nn)od(l1+1:)=offset(l2+1:)
        del=fit_mnb(start,yd,xd,n-1,mi,od);deallocate(od)
      else
        del=fit_mnb(start,yd,xd,n-1,mi)
      end if
      dpar=del%par-base%par;out%cook_distance(i)=dot_product(dpar,matmul(infoinv,dpar))
      out%likelihood_displacement(i)=2.0_dp*(base%loglik-del%loglik);deallocate(yd,xd)
    end do
  end function global_influence_mnb

  function local_influence_mnb(start,y,x,n,mi,scheme,covariate,offset) result(out)
    real(dp),intent(in)::start(:),y(:),x(:,:)
    integer,intent(in)::n,mi,scheme
    integer,intent(in),optional::covariate
    real(dp),intent(in),optional::offset(:)
    type(mnb_local_result)::out
    type(mnb_fit_result)::fit
    real(dp),allocatable::z(:),hh(:,:),delta(:,:),hinv(:,:),ff(:,:),eval(:),evec(:,:)
    integer::q,nw,jc,jc2
    logical::ok
    fit=fit_mnb(start,y,x,n,mi,offset);q=size(fit%par)
    select case(scheme)
    case(local_weight,local_dispersion,local_covariate);nw=n
    case(local_weight_obs);nw=n*mi
    case default;error stop 'local_influence_mnb: unknown scheme'
    end select
    jc=0;if(present(covariate))jc=covariate
    allocate(z(q+nw),hh(q+nw,q+nw),delta(q,nw),hinv(q,q),ff(nw,nw),eval(nw),evec(nw,nw))
    z(1:q)=fit%par
    if(scheme==local_covariate)then;z(q+1:)=0.0_dp;else;z(q+1:)=1.0_dp;end if
    call numerical_hessian(perturbed,z,hh);delta=hh(1:q,q+1:q+nw)
    call invert_matrix(fit%hessian,hinv,ok);if(.not.ok)error stop 'local_influence_mnb: singular Hessian'
    ff=matmul(transpose(delta),matmul(hinv,delta));call symmetric_eigen_jacobi(ff,eval,evec)
    allocate(out%direction(nw),out%total_curvature(nw));out%total_curvature=[(ff(jc2,jc2),jc2=1,nw)]
    ! R eigen() sorts decreasing; upstream selects the last eigenvector/value.
    out%selected_eigenvalue=eval(nw);out%direction=evec(:,nw)
  contains
    real(dp) function perturbed(zz) result(ll)
      real(dp),intent(in)::zz(:)
      real(dp),allocatable::eta(:),mu(:),ys(:),ms(:),xnew(:,:)
      real(dp)::phi,wi,sdx
      integer::i,k,nn,p,ii
      nn=n*mi;p=size(x,2);phi=zz(1);if(phi<=0.0_dp)then;ll=-huge(1.0_dp)/8.0_dp;return;end if
      allocate(eta(nn),mu(nn),ys(n),ms(n));eta=matmul(x,zz(2:q));if(present(offset))eta=eta+offset
      if(maxval(eta)>700.0_dp)then;ll=-huge(1.0_dp)/8.0_dp;return;end if
      mu=exp(eta);call mnb_cluster_sums(y,mu,n,mi,ys,ms);ll=0.0_dp
      select case(scheme)
      case(local_weight)
        do i=1,n
          wi=zz(q+i);ll=ll+wi*(log_gamma(phi+ys(i))-log_gamma(phi)+phi*log(phi)&
            -(phi+ys(i))*log(phi+ms(i)))
        end do
        do k=1,nn
          ii=mod(k-1,n)+1;wi=zz(q+ii);ll=ll+wi*(-log_gamma(y(k)+1.0_dp)+y(k)*eta(k))
        end do
      case(local_weight_obs)
        do k=1,nn
          i=(k-1)/mi+1;wi=zz(q+k)
          ll=ll+wi*((log_gamma(phi+ys(i))-log_gamma(phi)+phi*log(phi)-phi*log(phi+ms(i)))/real(mi,dp)&
            -log_gamma(y(k)+1.0_dp)+y(k)*eta(k)-y(k)*log(phi+ms(i)))
        end do
      case(local_covariate)
        if(jc<1 .or. jc>p)error stop 'local_influence_mnb: covariate index required'
        sdx=sample_sd_unique(x(:,jc));allocate(xnew(nn,p));xnew=x
        do k=1,nn;xnew(k,jc)=x(k,jc)+zz(q+mod(k-1,n)+1)*sdx;end do
        ! Upstream cova.pertu constructs X.new but then uses X in eta; preserve this exactly.
        ll=mnb_loglik(zz(1:q),y,x,n,mi,offset)
      case(local_dispersion)
        do i=1,n
          wi=zz(q+i)
          if(wi*phi<=0.0_dp)then;ll=-huge(1.0_dp)/8.0_dp;return;end if
          ll=ll+log_gamma(wi*phi+ys(i))-log_gamma(wi*phi)+wi*phi*log(wi*phi)&
            -(wi*phi+ys(i))*log(wi*phi+ms(i))
        end do
        ll=ll-sum(log_gamma(y+1.0_dp))+dot_product(y,eta)
      end select
    end function perturbed
  end function local_influence_mnb
end module mnb_influence
