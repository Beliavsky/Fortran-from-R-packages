! SPDX-License-Identifier: GPL-2.0-only
module ks_deconv
  use ks_kinds, only: dp
  use ks_normal, only: mvn_pdf
  use ks_kde, only: kde_model, fit_kde, kde_pdf
  implicit none
  private
  public :: deconv_weights, deconv_pdf, reg_ucv_value, reg_ucv
contains
  subroutine project_simplex(v,z)
    real(dp),intent(in)::v(:)
    real(dp),intent(out)::z(size(v))
    real(dp),allocatable::u(:)
    real(dp)::cssv,theta,tmp
    integer::i,j,rho,n
    n=size(v);allocate(u(n));u=v
    do i=2,n
      tmp=u(i);j=i-1
      do while(j>=1)
        if(u(j)>=tmp) exit
        u(j+1)=u(j);j=j-1
      end do
      u(j+1)=tmp
    end do
    cssv=0.0_dp;rho=1
    do i=1,n
      cssv=cssv+u(i)
      if(u(i)+(1.0_dp-cssv)/real(i,dp)>0.0_dp) rho=i
    end do
    theta=(sum(u(1:rho))-1.0_dp)/real(rho,dp)
    z=max(v-theta,0.0_dp)
    if(sum(z)>0.0_dp) z=z/sum(z)
  end subroutine

  real(dp) function spectral_upper(Q) result(lam)
    real(dp),intent(in)::Q(:,:)
    real(dp)::v(size(Q,1)),y(size(Q,1)),nv
    integer::it,n
    n=size(Q,1);v=1.0_dp/sqrt(real(n,dp));lam=1.0_dp
    do it=1,40
      y=matmul(Q,v);nv=sqrt(dot_product(y,y))
      if(nv<=tiny(1.0_dp)) then;lam=1.0_dp;return;end if
      v=y/nv;lam=max(dot_product(v,matmul(Q,v)),tiny(1.0_dp))
    end do
  end function

  subroutine deconv_weights(x,H,Sigma,reg,w,info,max_iter,tol)
    real(dp),intent(in)::x(:,:),H(:,:),Sigma(:,:),reg
    real(dp),allocatable,intent(out)::w(:)
    integer,intent(out),optional::info
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    real(dp),allocatable::Q(:,:),b(:),wn(:),grad(:)
    real(dp)::covq(size(H,1),size(H,2)),covb(size(H,1),size(H,2)),zero(size(H,1)),delta(size(H,1))
    real(dp)::step,eps,diff
    integer::n,d,i,j,it,nit
    n=size(x,1);d=size(x,2)
    if(n<=0.or.size(H,1)/=d.or.size(H,2)/=d.or.any(shape(Sigma)/=[d,d])) error stop 'deconv_weights: shape'
    allocate(Q(n,n),b(n),w(n),wn(n),grad(n));zero=0.0_dp
    covq=2.0_dp*H+2.0_dp*Sigma;covb=2.0_dp*H+Sigma
    do i=1,n
      do j=1,n
        delta=x(j,:)-x(i,:);Q(i,j)=mvn_pdf(delta,zero,covq)
      end do
      b(i)=0.0_dp
      do j=1,n
        delta=x(j,:)-x(i,:);b(i)=b(i)+mvn_pdf(delta,zero,covb)
      end do
      b(i)=b(i)/real(n,dp)
      Q(i,i)=Q(i,i)+reg/real(n,dp)
    end do
    w=1.0_dp/real(n,dp);step=1.0_dp/spectral_upper(Q)
    nit=10000;if(present(max_iter))nit=max_iter;eps=1.0e-10_dp;if(present(tol))eps=tol
    do it=1,nit
      grad=matmul(Q,w)-b
      call project_simplex(w-step*grad,wn)
      diff=maxval(abs(wn-w));w=wn
      if(diff<=eps) exit
    end do
    w=w*real(n,dp)
    if(present(info)) then
      if(it>nit) then;info=1;else;info=0;end if
    end if
  end subroutine

  subroutine deconv_pdf(x,H,Sigma,reg,eval,f,info)
    real(dp),intent(in)::x(:,:),H(:,:),Sigma(:,:),reg,eval(:,:)
    real(dp),intent(out)::f(size(eval,1))
    integer,intent(out),optional::info
    real(dp),allocatable::w(:)
    type(kde_model)::model
    integer::ierr
    call deconv_weights(x,H,Sigma,reg,w,ierr)
    call fit_kde(x,model,H=H,weights=w,info=ierr)
    call kde_pdf(model,eval,f)
    if(present(info))info=ierr
  end subroutine

  real(dp) function reg_ucv_value(x,H,Sigma,reg,kfold) result(cv)
    real(dp),intent(in)::x(:,:),H(:,:),Sigma(:,:),reg
    integer,intent(in),optional::kfold
    integer::k,n,fold,start,finish,ntr,nva,i,j,p
    real(dp),allocatable::train(:,:),valid(:,:),w(:),Ht(:,:),vals(:)
    type(kde_model)::model
    n=size(x,1);k=5;if(present(kfold))k=max(2,min(kfold,n))
    cv=0.0_dp
    do fold=1,k
      start=1+(fold-1)*n/k;finish=fold*n/k
      nva=finish-start+1;ntr=n-nva
      allocate(train(ntr,size(x,2)),valid(nva,size(x,2)),Ht(size(H,1),size(H,2)),vals(nva))
      i=0;p=0
      do j=1,n
        if(j>=start.and.j<=finish)then;p=p+1;valid(p,:)=x(j,:)
        else;i=i+1;train(i,:)=x(j,:)
        end if
      end do
      call deconv_weights(train,H,Sigma,reg,w)
      Ht=H+Sigma
      call fit_kde(train,model,H=Ht,weights=w)
      call kde_pdf(model,valid,vals);cv=cv+sum(vals)
      deallocate(train,valid,Ht,vals,w)
    end do
  end function

  real(dp) function reg_ucv(x,H,Sigma,kfold) result(reg)
    real(dp),intent(in)::x(:,:),H(:,:),Sigma(:,:)
    integer,intent(in),optional::kfold
    real(dp)::a,b,c,d,fc,fd,phi
    integer::it
    phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp;a=-10.0_dp;b=4.0_dp
    c=b-phi*(b-a);d=a+phi*(b-a)
    fc=-reg_ucv_value(x,H,Sigma,exp(2.0_dp*c),kfold);fd=-reg_ucv_value(x,H,Sigma,exp(2.0_dp*d),kfold)
    do it=1,60
      if(fc<fd)then;b=d;d=c;fd=fc;c=b-phi*(b-a);fc=-reg_ucv_value(x,H,Sigma,exp(2.0_dp*c),kfold)
      else;a=c;c=d;fc=fd;d=a+phi*(b-a);fd=-reg_ucv_value(x,H,Sigma,exp(2.0_dp*d),kfold)
      end if
    end do
    reg=exp(2.0_dp*0.5_dp*(a+b))
  end function
end module ks_deconv
