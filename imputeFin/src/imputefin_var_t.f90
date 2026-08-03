! SPDX-License-Identifier: GPL-3.0-only
module imputefin_var_t
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use imputefin_kinds, only : dp
  use imputefin_types, only : var_t_options, var_t_result, impute_ok, impute_invalid_input, &
       impute_insufficient_data, impute_singular, impute_not_converged
  use imputefin_linalg, only : solve_spd, inverse_matrix, symmetrize
  use imputefin_math, only : log_gamma_dp
  implicit none
  private
  public :: fit_var_t
contains
  pure logical function close_rel_array(a,b,tol)
    real(dp), intent(in) :: a(:),b(:),tol
    close_rel_array=all(abs(a-b)<=tol*max(1.0_dp,0.5_dp*(abs(a)+abs(b))))
  end function close_rel_array

  subroutine set_error(res,status,message)
    type(var_t_result),intent(inout)::res
    integer,intent(in)::status
    character(*),intent(in)::message
    res%status=status;res%message=message;res%converged=.false.
  end subroutine set_error

  subroutine fill_column_means(y)
    real(dp),intent(inout)::y(:,:)
    integer::j,i,nobs
    real(dp)::mu
    do j=1,size(y,2)
      mu=0.0_dp;nobs=0
      do i=1,size(y,1)
        if(.not.ieee_is_nan(y(i,j)))then;mu=mu+y(i,j);nobs=nobs+1;end if
      end do
      if(nobs>0)mu=mu/real(nobs,dp)
      do i=1,size(y,1);if(ieee_is_nan(y(i,j)))y(i,j)=mu;end do
    end do
  end subroutine fill_column_means

  subroutine build_design(y,p,x,z,rows,complete_only,template,status)
    real(dp),intent(in)::y(:,:),template(:,:)
    integer,intent(in)::p
    real(dp),allocatable,intent(out)::x(:,:),z(:,:)
    integer,allocatable,intent(out)::rows(:)
    logical,intent(in)::complete_only
    integer,intent(out)::status
    logical::ok
    integer::t,j,lag,k,m,n,q
    n=size(y,1);q=size(y,2);m=0
    do t=p+1,n
      ok=.true.
      if(complete_only)then
        do j=1,q
          if(ieee_is_nan(template(t,j)))ok=.false.
          do lag=1,p;if(ieee_is_nan(template(t-lag,j)))ok=.false.;end do
        end do
      end if
      if(ok)m=m+1
    end do
    if(m<max(5,1+p*q))then;status=1;allocate(x(0,0),z(0,0),rows(0));return;end if
    allocate(x(m,1+p*q),z(m,q),rows(m));k=0
    do t=p+1,n
      ok=.true.
      if(complete_only)then
        do j=1,q
          if(ieee_is_nan(template(t,j)))ok=.false.
          do lag=1,p;if(ieee_is_nan(template(t-lag,j)))ok=.false.;end do
        end do
      end if
      if(.not.ok)cycle
      k=k+1;rows(k)=t;x(k,1)=1.0_dp;z(k,:)=y(t,:)
      do lag=1,p
        x(k,2+(lag-1)*q:1+lag*q)=y(t-lag,:)
      end do
    end do
    status=0
  end subroutine build_design

  subroutine weighted_var_fit(x,z,w,b,scatter,status)
    real(dp),intent(in)::x(:,:),z(:,:),w(:)
    real(dp),intent(out)::b(:,:),scatter(:,:)
    integer,intent(out)::status
    real(dp),allocatable::xtwx(:,:),xtwz(:,:),sol(:,:),e(:,:)
    integer::i,m,q,info
    m=size(x,1);q=size(z,2)
    allocate(xtwx(size(x,2),size(x,2)),xtwz(size(x,2),q),sol(size(x,2),q),e(m,q))
    xtwx=0.0_dp;xtwz=0.0_dp
    do i=1,m
      xtwx=xtwx+w(i)*spread(x(i,:),2,size(x,2))*spread(x(i,:),1,size(x,2))
      xtwz=xtwz+w(i)*spread(x(i,:),2,q)*spread(z(i,:),1,size(x,2))
    end do
    do i=1,size(xtwx,1);xtwx(i,i)=xtwx(i,i)+1.0e-10_dp;end do
    call solve_spd(xtwx,xtwz,sol,info)
    if(info/=0)then;status=1;return;end if
    b=transpose(sol);e=z-matmul(x,transpose(b));scatter=0.0_dp
    do i=1,m;scatter=scatter+w(i)*spread(e(i,:),2,q)*spread(e(i,:),1,q);end do
    scatter=scatter/real(m,dp);call symmetrize(scatter)
    do i=1,q;scatter(i,i)=max(scatter(i,i),1.0e-10_dp);end do
    status=0
  end subroutine weighted_var_fit

  function nu_obj(nu,delta,q) result(v)
    real(dp),intent(in)::nu,delta(:)
    integer,intent(in)::q
    real(dp)::v
    if(nu<=2.000001_dp)then;v=huge(1.0_dp);return;end if
    v=sum(0.5_dp*(nu+real(q,dp))*log(1.0_dp+delta/nu))+real(size(delta),dp)*(&
      log_gamma_dp(0.5_dp*nu)-log_gamma_dp(0.5_dp*(nu+real(q,dp)))+0.5_dp*real(q,dp)*log(nu))
  end function nu_obj

  subroutine optimize_nu(delta,q,nu)
    real(dp),intent(in)::delta(:)
    integer,intent(in)::q
    real(dp),intent(out)::nu
    real(dp)::a,b,x1,x2,f1,f2,gr
    integer::it
    a=2.0001_dp;b=100.0_dp;gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    x1=b-gr*(b-a);x2=a+gr*(b-a);f1=nu_obj(x1,delta,q);f2=nu_obj(x2,delta,q)
    do it=1,150
      if(abs(b-a)<1.0e-7_dp*(1.0_dp+abs(x1)+abs(x2)))exit
      if(f1>f2)then;a=x1;x1=x2;f1=f2;x2=a+gr*(b-a);f2=nu_obj(x2,delta,q)
      else;b=x2;x2=x1;f2=f1;x1=b-gr*(b-a);f1=nu_obj(x1,delta,q);end if
    end do
    nu=merge(x1,x2,f1<=f2)
  end subroutine optimize_nu

  subroutine unpack_b(b,p,phi0,phi)
    real(dp),intent(in)::b(:,:)
    integer,intent(in)::p
    real(dp),intent(out)::phi0(:),phi(:,:,:)
    integer::lag,q
    q=size(b,1);phi0=b(:,1)
    do lag=1,p;phi(:,:,lag)=b(:,2+(lag-1)*q:1+lag*q);end do
  end subroutine unpack_b

  subroutine conditional_fill_rows(y,template,p,b,scatter)
    real(dp),intent(inout)::y(:,:)
    real(dp),intent(in)::template(:,:),b(:,:),scatter(:,:)
    integer,intent(in)::p
    integer::t,j,lag,q,nm,no,info
    logical,allocatable::miss(:),obs(:)
    integer,allocatable::im(:),io(:)
    real(dp),allocatable::xp(:),mu(:),soo(:,:),smo(:,:),rhs(:,:),sol(:,:)
    q=size(y,2);allocate(xp(1+p*q),mu(q),miss(q),obs(q))
    do t=p+1,size(y,1)
      xp(1)=1.0_dp
      do lag=1,p;xp(2+(lag-1)*q:1+lag*q)=y(t-lag,:);end do
      mu=matmul(b,xp)
      do j=1,q;miss(j)=ieee_is_nan(template(t,j));obs(j)=.not.miss(j);end do
      if(.not.any(miss))cycle
      call mask_indices(miss,im);call mask_indices(obs,io);nm=size(im);no=size(io)
      if(no==0)then;y(t,im)=mu(im)
      else
        soo=scatter(io,io);smo=scatter(im,io);allocate(rhs(no,1),sol(no,1));rhs(:,1)=template(t,io)-mu(io)
        call solve_spd(soo,rhs,sol,info)
        if(info==0)y(t,im)=mu(im)+matmul(smo,sol(:,1))
        deallocate(rhs,sol)
      end if
      if(allocated(im))deallocate(im);if(allocated(io))deallocate(io)
      if(allocated(soo))deallocate(soo);if(allocated(smo))deallocate(smo)
    end do
  end subroutine conditional_fill_rows

  subroutine mask_indices(mask,idx)
    logical,intent(in)::mask(:)
    integer,allocatable,intent(out)::idx(:)
    integer::i,k
    allocate(idx(count(mask)));k=0
    do i=1,size(mask);if(mask(i))then;k=k+1;idx(k)=i;end if;end do
  end subroutine mask_indices

  subroutine fit_var_t(y,res,options)
    real(dp),intent(in)::y(:,:)
    type(var_t_result),intent(out)::res
    type(var_t_options),intent(in),optional::options
    type(var_t_options)::opt
    real(dp),allocatable::work(:,:),x(:,:),z(:,:),w(:),b(:,:),bold(:,:),scatter(:,:),sold(:,:),&
         invs(:,:),e(:,:),delta(:),rhs(:,:)
    integer,allocatable::rows(:)
    integer::p,q,m,k,k_used,i,j,info,st
    logical :: had_missing
    real(dp)::nu,nuold,nunew
    opt=var_t_options();if(present(options))opt=options
    res=var_t_result();p=opt%p;q=size(y,2)
    if(p<1.or.size(y,1)<=p+4.or.q<1.or.opt%maxiter<1.or.opt%tol<=0.0_dp)then
      call set_error(res,impute_invalid_input,'invalid dimensions or options');return
    end if
    had_missing=.false.
    do j=1,size(y,2)
      do i=1,size(y,1)
        if(ieee_is_nan(y(i,j)))had_missing=.true.
      end do
    end do
    work=y
    call fill_column_means(work)
    call build_design(work,p,x,z,rows,opt%omit_missing,y,st)
    if(st/=0)then;call set_error(res,impute_insufficient_data,'too few usable VAR rows');return;end if
    m=size(x,1);allocate(w(m),b(q,size(x,2)),bold(q,size(x,2)),scatter(q,q),sold(q,q),delta(m),e(m,q),invs(q,q),rhs(q,1))
    w=1.0_dp;call weighted_var_fit(x,z,w,b,scatter,st)
    if(st/=0)then;call set_error(res,impute_singular,'initial VAR regression failed');return;end if
    nu=6.0_dp;res%converged=.false.
    do k=1,opt%maxiter
      bold=b;sold=scatter;nuold=nu
      if(.not.opt%omit_missing .and. had_missing)then
        call conditional_fill_rows(work,y,p,b,scatter)
        call build_design(work,p,x,z,rows,.false.,y,st);m=size(x,1)
        if(size(w)/=m)then;deallocate(w,delta,e);allocate(w(m),delta(m),e(m,q));end if
      end if
      call inverse_matrix(scatter,invs,info)
      if(info/=0)then
        call set_error(res,impute_singular,'scatter inverse failed')
        return
      end if
      e=z-matmul(x,transpose(b))
      do i=1,m
        rhs(:,1)=e(i,:);delta(i)=dot_product(rhs(:,1),matmul(invs,rhs(:,1)))
      end do
      w=(nu+real(q,dp))/(nu+delta)
      call weighted_var_fit(x,z,w,b,scatter,st)
      if(st/=0)then
        call set_error(res,impute_singular,'weighted VAR regression failed')
        return
      end if
      call optimize_nu(delta,q,nunew);nu=nunew
      if(close_rel_array(reshape(b,[size(b)]),reshape(bold,[size(bold)]),opt%tol).and.&
         close_rel_array(reshape(scatter,[q*q]),reshape(sold,[q*q]),opt%tol).and.&
         abs(nu-nuold)<=sqrt(opt%tol)*max(1.0_dp,0.5_dp*(nu+nuold)))then;res%converged=.true.;exit;end if
    end do
    allocate(res%phi0(q),res%phi(q,q,p),res%scatter(q,q),res%completed(size(work,1),q))
    k_used=min(k,opt%maxiter)
    call unpack_b(b,p,res%phi0,res%phi)
    res%scatter=scatter
    res%completed=work
    res%nu=nu
    res%iterations=k_used
    res%n_used=m
    res%status=merge(impute_ok,impute_not_converged,res%converged)
    if (res%converged) then
      res%message='ok'
    else
      res%message='maximum iterations reached'
    end if
  end subroutine fit_var_t
end module imputefin_var_t
