module l1pack_multivariate
  use l1pack_base, only: dp, pi, chol_lower, mahalanobis_one, weighted_center_scatter, &
    weighted_center, as78_median_center, bessel_k_ratio_prev, logdet_chol
  use l1pack_distributions, only: log_dmlaplace, rmlaplace
  implicit none
  private

  type, public :: laplace_fit_result
    real(dp), allocatable :: center(:), scatter(:,:), distances(:), weights(:)
    real(dp) :: loglik=0.0_dp
    integer :: n=0,p=0,iterations=0
    logical :: converged=.false.
  end type laplace_fit_result

  type, public :: spatial_median_result
    real(dp), allocatable :: median(:), scatter(:,:), distances(:), weights(:)
    real(dp) :: loglik=0.0_dp
    integer :: n=0,p=0,iterations=0,inner_iterations=0
    logical :: converged=.false.
  end type spatial_median_result

  type, public :: laplace_envelope_result
    real(dp), allocatable :: transformed(:), lower(:), upper(:)
  end type laplace_envelope_result

  public :: laplace_fit, laplace_fit_equal, spatial_median_fit
  public :: wh_laplace, envelope_laplace

contains

  subroutine initial_center_scatter(x,center,scatter)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::center(:),scatter(:,:)
    real(dp)::w(size(x,1))
    w=1.0_dp
    call weighted_center_scatter(x,w,center,scatter)
    call ensure_pd(scatter)
  end subroutine initial_center_scatter

  subroutine ensure_pd(a)
    real(dp),intent(inout)::a(:,:)
    real(dp)::l(size(a,1),size(a,2)),jitter,tr
    integer::info,i,k
    a=0.5_dp*(a+transpose(a))
    call chol_lower(a,l,info)
    if(info==0)return
    tr=sum([(a(i,i),i=1,size(a,1))])/real(size(a,1),dp)
    jitter=max(1.0e-10_dp,1.0e-8_dp*max(abs(tr),1.0_dp))
    do k=1,12
      do i=1,size(a,1);a(i,i)=a(i,i)+jitter;end do
      call chol_lower(a,l,info)
      if(info==0)return
      jitter=jitter*10.0_dp
    end do
  end subroutine ensure_pd

  subroutine distances_from_scatter(x,center,scatter,dist)
    real(dp),intent(in)::x(:,:),center(:),scatter(:,:)
    real(dp),intent(out)::dist(:)
    integer::i
    do i=1,size(x,1)
      dist(i)=mahalanobis_one(x(i,:),center,scatter)
    end do
  end subroutine distances_from_scatter

  real(dp) function multiv_laplace_loglik(x,center,scatter) result(v)
    real(dp),intent(in)::x(:,:),center(:),scatter(:,:)
    integer::i
    v=0.0_dp
    do i=1,size(x,1);v=v+log_dmlaplace(x(i,:),center,scatter);end do
  end function multiv_laplace_loglik

  real(dp) function do_weight(p,d2) result(w)
    integer,intent(in)::p
    real(dp),intent(in)::d2
    real(dp)::d,x,ratio,nu
    d=sqrt(max(d2,0.0_dp))
    if(d<=1.0e-12_dp) then
      w=1.0e12_dp
      return
    end if
    x=0.5_dp*d;nu=0.5_dp*real(p,dp)
    ratio=bessel_k_ratio_prev(nu,x)
    w=0.5_dp*ratio/d
    w=min(w,1.0e12_dp)
  end function do_weight

  subroutine laplace_fit(x,res,tol,maxiter)
    real(dp),intent(in)::x(:,:)
    type(laplace_fit_result),intent(out)::res
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer::n,p,i,iter,mi
    real(dp)::tt,oldll,newll,conv
    real(dp),allocatable::center(:),scatter(:,:),dist(:),weights(:)

    n=size(x,1);p=size(x,2)
    if(n<2.or.p<1)error stop 'laplace_fit: invalid data dimensions'
    tt=1.0e-6_dp;if(present(tol))tt=tol
    mi=200;if(present(maxiter))mi=maxiter
    allocate(center(p),scatter(p,p),dist(n),weights(n))
    call initial_center_scatter(x,center,scatter)
    call distances_from_scatter(x,center,scatter,dist)
    oldll=multiv_laplace_loglik(x,center,scatter)

    do iter=1,mi
      do i=1,n;weights(i)=do_weight(p,dist(i));end do
      call weighted_center_scatter(x,weights,center,scatter)
      call ensure_pd(scatter)
      call distances_from_scatter(x,center,scatter,dist)
      newll=multiv_laplace_loglik(x,center,scatter)
      conv=abs((newll-oldll)/(abs(newll)+1.0e-50_dp))
      if(conv<tt)exit
      oldll=newll
    end do
    res%n=n;res%p=p;res%iterations=min(iter,mi);res%converged=(iter<=mi)
    allocate(res%center(p),res%scatter(p,p),res%distances(n),res%weights(n))
    res%center=center;res%scatter=scatter;res%distances=sqrt(max(dist,0.0_dp));res%weights=weights
    res%loglik=multiv_laplace_loglik(x,center,scatter)
  end subroutine laplace_fit

  subroutine laplace_fit_equal(x,start,res,tol,maxiter)
    real(dp),intent(in)::x(:,:)
    type(laplace_fit_result),intent(in)::start
    type(laplace_fit_result),intent(out)::res
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer::n,p,i,j,iter,mi
    real(dp)::tt,oldll,newll,conv,lambda
    real(dp),allocatable::center(:),meanvec(:),scatter(:,:),dist(:),weights(:),l(:,:),ones(:),z(:),tmp(:)
    integer::info

    n=size(x,1);p=size(x,2);tt=1.0e-6_dp;if(present(tol))tt=tol;mi=200;if(present(maxiter))mi=maxiter
    allocate(center(p),meanvec(p),scatter(p,p),dist(n),weights(n),l(p,p),ones(p),z(p),tmp(p))
    center=start%center;scatter=start%scatter;lambda=sum(center)/real(p,dp);meanvec=lambda
    call distances_from_scatter(x,meanvec,scatter,dist);oldll=multiv_laplace_loglik(x,meanvec,scatter)
    do iter=1,mi
      do i=1,n;weights(i)=do_weight(p,dist(i));end do
      call weighted_center(x,weights,center)
      call chol_lower(scatter,l,info);if(info/=0)exit
      ones=1.0_dp;z=center
      tmp=matmul(l,ones);z=matmul(l,z)
      lambda=dot_product(tmp,z)/dot_product(tmp,tmp)
      meanvec=lambda
      scatter=0.0_dp
      do i=1,n
        z=x(i,:)-lambda
        do j=1,p;scatter(j,:)=scatter(j,:)+weights(i)*z(j)*z;end do
      end do
      scatter=scatter/real(n,dp);call ensure_pd(scatter)
      call distances_from_scatter(x,meanvec,scatter,dist)
      newll=multiv_laplace_loglik(x,meanvec,scatter)
      conv=abs((newll-oldll)/(abs(newll)+1.0e-50_dp))
      if(conv<tt)exit
      oldll=newll
    end do
    res%n=n;res%p=p;res%iterations=min(iter,mi);res%converged=(iter<=mi)
    allocate(res%center(p),res%scatter(p,p),res%distances(n),res%weights(n))
    res%center=meanvec;res%scatter=scatter;res%distances=sqrt(max(dist,0.0_dp));res%weights=weights
    res%loglik=multiv_laplace_loglik(x,meanvec,scatter)
  end subroutine laplace_fit_equal

  subroutine forward_solve(l,b,x)
    real(dp),intent(in)::l(:,:),b(:)
    real(dp),intent(out)::x(:)
    integer::i
    x=0.0_dp
    do i=1,size(b)
      x(i)=(b(i)-dot_product(l(i,1:i-1),x(1:i-1)))/l(i,i)
    end do
  end subroutine forward_solve

  subroutine spatial_median_fit(x,res,tol,maxiter)
    real(dp),intent(in)::x(:,:)
    type(spatial_median_result),intent(out)::res
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer::n,p,i,j,iter,mi,inner,totinner,info
    real(dp)::tt,oldll,newll,conv,d
    real(dp),allocatable::med(:),scatter(:,:),dist(:),weights(:),l(:,:),z(:,:),tmp(:)

    n=size(x,1);p=size(x,2);tt=1.0e-6_dp;if(present(tol))tt=tol;mi=200;if(present(maxiter))mi=maxiter
    allocate(med(p),scatter(p,p),dist(n),weights(n),l(p,p),z(n,p),tmp(p))
    call initial_center_scatter(x,med,scatter)
    call distances_from_scatter(x,med,scatter,dist)
    oldll=kotz_loglik(x,med,scatter);totinner=0
    do iter=1,mi
      call chol_lower(scatter,l,info);if(info/=0)exit
      do i=1,n
        call forward_solve(l,x(i,:),tmp);z(i,:)=tmp
      end do
      call as78_median_center(z,tmp,inner);totinner=totinner+abs(inner)
      med=matmul(l,tmp)
      call distances_from_scatter(x,med,scatter,dist)
      do i=1,n
        d=sqrt(max(dist(i),1.0e-24_dp));weights(i)=1.0_dp/d
      end do
      scatter=0.0_dp
      do i=1,n
        tmp=x(i,:)-med
        do j=1,p;scatter(j,:)=scatter(j,:)+weights(i)*tmp(j)*tmp;end do
      end do
      scatter=scatter/real(n,dp);call ensure_pd(scatter)
      call distances_from_scatter(x,med,scatter,dist)
      newll=kotz_loglik(x,med,scatter)
      conv=abs((newll-oldll)/(abs(newll)+1.0e-50_dp))
      if(conv<tt)exit
      oldll=newll
    end do
    res%n=n;res%p=p;res%iterations=min(iter,mi);res%inner_iterations=totinner;res%converged=(iter<=mi)
    allocate(res%median(p),res%scatter(p,p),res%distances(n),res%weights(n))
    res%median=med;res%scatter=scatter;res%distances=sqrt(max(dist,0.0_dp));res%weights=weights
    res%loglik=kotz_loglik(x,med,scatter)
  end subroutine spatial_median_fit

  real(dp) function kotz_loglik(x,med,scatter) result(v)
    real(dp),intent(in)::x(:,:),med(:),scatter(:,:)
    integer::n,p,i,info
    real(dp)::ld,acc
    n=size(x,1);p=size(x,2);ld=logdet_chol(scatter,info)
    if(info/=0)then;v=-huge(1.0_dp);return;end if
    acc=0.0_dp
    do i=1,n;acc=acc+sqrt(max(mahalanobis_one(x(i,:),med,scatter),0.0_dp));end do
    v=real(n,dp)*(-real(p,dp)*log(2.0_dp)-0.5_dp*real(p-1,dp)*log(pi) &
      -log_gamma(0.5_dp*real(p+1,dp))-ld)-acc
  end function kotz_loglik

  pure function wh_laplace(distances,p) result(z)
    real(dp),intent(in)::distances(:)
    integer,intent(in)::p
    real(dp)::z(size(distances)),f,meanv,sd
    integer::i
    meanv=1.0_dp-1.0_dp/(9.0_dp*real(p,dp));sd=1.0_dp/sqrt(9.0_dp*real(p,dp))
    do i=1,size(distances)
      f=distances(i)/(2.0_dp*real(p,dp))
      z(i)=(max(f,0.0_dp)**(1.0_dp/3.0_dp)-meanv)/sd
    end do
  end function wh_laplace

  subroutine envelope_laplace(x,fit,reps,conf,out)
    real(dp),intent(in)::x(:,:)
    type(laplace_fit_result),intent(in)::fit
    integer,intent(in),optional::reps
    real(dp),intent(in),optional::conf
    type(laplace_envelope_result),intent(out)::out
    integer::n,p,rp,b,i,lo,hi
    real(dp)::cf
    real(dp),allocatable::sim(:,:),z(:),elims(:,:),row(:)
    type(laplace_fit_result)::sf
    n=size(x,1);p=size(x,2);rp=50;if(present(reps))rp=reps;cf=0.95_dp;if(present(conf))cf=conf
    allocate(out%transformed(n),out%lower(n),out%upper(n),sim(n,p),z(n),elims(n,rp),row(rp))
    out%transformed=wh_laplace(fit%distances,p)
    do b=1,rp
      call rmlaplace(n,fit%center,fit%scatter,sim)
      call laplace_fit(sim,sf)
      z=wh_laplace(sf%distances,p)
      call sort_inplace(z);elims(:,b)=z
    end do
    lo=max(1,int((1.0_dp-cf)*0.5_dp*real(rp,dp))+1)
    hi=min(rp,int((1.0_dp-(1.0_dp-cf)*0.5_dp)*real(rp,dp)))
    do i=1,n
      row=elims(i,:);call sort_inplace(row);out%lower(i)=row(lo);out%upper(i)=row(max(lo,hi))
    end do
  end subroutine envelope_laplace

  subroutine sort_inplace(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::t
    do i=2,size(x)
      t=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=t)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=t
    end do
  end subroutine sort_inplace

end module l1pack_multivariate
