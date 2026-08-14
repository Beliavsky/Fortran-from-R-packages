! SPDX-License-Identifier: GPL-2.0-only
module ks_feature
  use ks_kinds, only: dp
  use ks_kde, only: kde_model, kde_pdf, kdde_eval
  use ks_linalg, only: symmetric_eigen, spd_inverse
  use ks_normal, only: mvn_derivative_tensor
  use ks_special, only: chi_square_cdf
  implicit none
  private
  public :: kcurv_eval, kfs_eval
contains
  subroutine decode_hessian(v,d,a)
    real(dp),intent(in)::v(:)
    integer,intent(in)::d
    real(dp),intent(out)::a(d,d)
    integer::i,j
    if(size(v)/=d*d) error stop 'decode_hessian: shape'
    do j=1,d;do i=1,d;a(i,j)=v(i+(j-1)*d);end do;end do
    a=0.5_dp*(a+transpose(a))
  end subroutine

  subroutine kcurv_eval(model,points,curvature,local_mode)
    type(kde_model),intent(in)::model
    real(dp),intent(in)::points(:,:)
    real(dp),intent(out)::curvature(size(points,1))
    logical,intent(out),optional::local_mode(size(points,1))
    real(dp),allocatable::d2(:,:)
    real(dp)::hess(model%d,model%d),eig(model%d),vecs(model%d,model%d)
    integer::i,info
    logical::lm
    call kdde_eval(model,points,2,d2)
    do i=1,size(points,1)
      call decode_hessian(d2(i,:),model%d,hess)
      call symmetric_eigen(hess,eig,vecs,info)
      lm=info==0 .and. all(eig<=0.0_dp)
      if(lm) then;curvature(i)=abs(product(eig));else;curvature(i)=0.0_dp;end if
      if(present(local_mode))local_mode(i)=lm
    end do
  end subroutine

  subroutine kfs_eval(model,points,signif_level,significant,wald,pvalue,local_mode)
    type(kde_model),intent(in)::model
    real(dp),intent(in)::points(:,:),signif_level
    logical,intent(out)::significant(size(points,1))
    real(dp),intent(out),optional::wald(size(points,1)),pvalue(size(points,1))
    logical,intent(out),optional::local_mode(size(points,1))
    integer::d,q,p,r,i,j,k,l,info,idx4,m,rankpos,ntest
    integer,allocatable::pi1(:),pi2(:),ord(:)
    real(dp),allocatable::d2(:,:),d4(:),C(:,:),Cinv(:,:),hv(:),dens(:),wv(:),pv(:)
    real(dp)::zero(model%d),cov2(model%d,model%d),scale,stat,thresh
    real(dp)::hess(model%d,model%d),eig(model%d),evec(model%d,model%d)
    logical,allocatable::lm(:),candidate(:)
    d=model%d;q=d*(d+1)/2
    allocate(pi1(q),pi2(q),C(q,q),Cinv(q,q),hv(q),dens(size(points,1)),wv(size(points,1)), &
             pv(size(points,1)),lm(size(points,1)),candidate(size(points,1)),ord(size(points,1)))
    p=0
    do j=1,d;do i=j,d;p=p+1;pi1(p)=i;pi2(p)=j;end do;end do
    zero=0.0_dp;cov2=2.0_dp*model%H
    call mvn_derivative_tensor(zero,zero,cov2,4,d4,info)
    if(info/=0) error stop 'kfs_eval: derivative covariance'
    do p=1,q
      do r=1,q
        i=pi1(p);j=pi2(p);k=pi1(r);l=pi2(r)
        idx4=1+(i-1)+d*(j-1)+d*d*(k-1)+d*d*d*(l-1)
        C(p,r)=d4(idx4)
      end do
    end do
    scale=sum(model%w*model%w)/real(model%n*model%n,dp)
    C=C*scale
    call spd_inverse(C,Cinv,info)
    if(info/=0) error stop 'kfs_eval: covariance singular'
    call kde_pdf(model,points,dens);call kdde_eval(model,points,2,d2)
    do m=1,size(points,1)
      call decode_hessian(d2(m,:),d,hess)
      call symmetric_eigen(hess,eig,evec,info)
      lm(m)=info==0 .and. all(eig<=0.0_dp)
      p=0
      do j=1,d;do i=j,d;p=p+1;hv(p)=hess(i,j);end do;end do
      if(dens(m)>tiny(1.0_dp))then;stat=dot_product(hv,matmul(Cinv,hv))/dens(m);else;stat=huge(1.0_dp);end if
      wv(m)=stat;pv(m)=1.0_dp-chi_square_cdf(stat,q)
    end do
    candidate=lm
    ! Hochberg step-up adjustment over candidate local modes.
    ord=[(m,m=1,size(points,1))]
    do i=2,size(ord)
      m=ord(i);j=i-1
      do while(j>=1)
        if(pv(ord(j))<=pv(m))exit
        ord(j+1)=ord(j);j=j-1
      end do
      ord(j+1)=m
    end do
    ntest=count(candidate);significant=.false.;rankpos=0
    do i=1,size(ord)
      m=ord(i);if(.not.candidate(m))cycle
      rankpos=rankpos+1;thresh=signif_level/real(ntest+1-rankpos,dp)
      if(pv(m)<=thresh)significant(m)=.true.
    end do
    if(present(wald))wald=wv;if(present(pvalue))pvalue=pv;if(present(local_mode))local_mode=lm
  end subroutine
end module ks_feature
