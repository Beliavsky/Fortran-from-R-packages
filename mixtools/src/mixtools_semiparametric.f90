! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_semiparametric
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mixtools_kinds, only : dp
  use mixtools_status
  use mixtools_types
  use mixtools_distributions, only : normal_pdf, normalize_logweights
  use mixtools_utilities, only : wkde, weighted_bandwidth, sort_real
  use mixtools_regression, only : regmix_em
  implicit none
  private
  public :: npem, npem_indrep, npem_indrepbw, npmsl, spem, mvnpem
  public :: spem_symloc, spem_symloc_n01, spregmix
contains
  subroutine make_grid(x,ngrid,grid)
    real(dp),intent(in)::x(:)
    integer,intent(in)::ngrid
    real(dp),allocatable,intent(out)::grid(:)
    real(dp)::lo,hi,pad
    integer::i
    lo=minval(x);hi=maxval(x);pad=max(0.05_dp*(hi-lo),1.0e-3_dp)
    allocate(grid(ngrid));do i=1,ngrid;grid(i)=lo-pad+(hi-lo+2.0_dp*pad)*real(i-1,dp)/real(max(1,ngrid-1),dp);end do
  end subroutine make_grid

  subroutine initialize_post(x,k,post)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::k
    real(dp),intent(out)::post(size(x,1),k)
    real(dp),allocatable::score(:),sorted(:)
    real(dp)::scale,center
    integer::i,j,n,idx
    n=size(x,1);allocate(score(n),sorted(n));score=sum(x,dim=2)/real(size(x,2),dp);sorted=score;call sort_real(sorted)
    scale=sqrt(sum((score-sum(score)/real(n,dp))**2)/real(max(1,n-1),dp));scale=max(scale,1.0e-3_dp)
    do j=1,k
      idx=max(1,min(n,nint((real(j,dp)-0.5_dp)*real(n,dp)/real(k,dp))));center=sorted(idx)
      post(:,j)=exp(-0.5_dp*((score-center)/scale)**2)
    end do
    do i=1,n;post(i,:)=post(i,:)/sum(post(i,:));end do
  end subroutine initialize_post

  subroutine product_kde_em(x,k,result,control,bandwidth,same_bandwidth,ngrid,symmetric)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::k
    type(semiparametric_result),intent(out)::result
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:)
    logical,intent(in),optional::same_bandwidth,symmetric
    integer,intent(in),optional::ngrid
    type(em_control)::ctl
    integer::n,p,m,i,j,q,iter,status
    logical::samebw,sym
    real(dp),allocatable::post(:,:),lambda(:),bw(:),history(:),lw(:),densobs(:,:,:),grid(:),densgrid(:,:),tmp(:)
    real(dp)::ll,newll,ln,diff,meanbw
    ctl=em_control();if(present(control))ctl=control;samebw=.true.;sym=.false.;m=200
    if(present(same_bandwidth))samebw=same_bandwidth;if(present(symmetric))sym=symmetric;if(present(ngrid))m=ngrid
    n=size(x,1);p=size(x,2)
    if(n<k.or.k<1)then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(post(n,k),lambda(k),bw(p),history(ctl%max_iterations+1),lw(k),densobs(n,p,k))
    call initialize_post(x,k,post);lambda=sum(post,dim=1)/real(n,dp)
    if(present(bandwidth))then
      if(size(bandwidth)==1)bw=bandwidth(1)
      if(size(bandwidth)==p)bw=bandwidth
    else
      do q=1,p;bw(q)=weighted_bandwidth(x(:,q),spread(1.0_dp,1,n));end do
    end if
    if(samebw)then;meanbw=sum(bw)/real(p,dp);bw=meanbw;end if
    call update_densities();call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      lambda=max(sum(post,dim=1)/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      call update_densities();call estep(newll);if(status/=0)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    call make_grid(reshape(x,[n*p]),m,grid);allocate(densgrid(m,k),tmp(m));densgrid=0.0_dp
    do j=1,k
      do q=1,p;call wkde(x(:,q),post(:,j),grid,bw(q),tmp,sym);densgrid(:,j)=densgrid(:,j)+tmp/real(p,dp);end do
    end do
    result%lambda=lambda;result%posterior=post;result%grid=grid;result%density=transpose(densgrid)
    result%bandwidth=bw;result%loglik=ll;result%iterations=min(iter,ctl%max_iterations)
    result%loglik_history=history(:result%iterations+1);result%converged=iter<=ctl%max_iterations.and.status==0
    result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine update_densities()
      real(dp),allocatable::vals(:)
      allocate(vals(n))
      do j=1,k;do q=1,p;call wkde(x(:,q),post(:,j),x(:,q),bw(q),vals,sym);densobs(:,q,j)=max(vals,tiny(1.0_dp));end do;end do
    end subroutine update_densities
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do i=1,n;do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+sum(log(densobs(i,:,j)));end do
        call normalize_logweights(lw,post(i,:),ln);if(.not.ieee_is_finite(ln))then;status=MIXTOOLS_NUMERICAL_ERROR;return;end if
        ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine product_kde_em

  subroutine npem(x,k,result,control,bandwidth,same_bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);logical,intent(in),optional::same_bandwidth;integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,same_bandwidth,ngrid,.false.)
  end subroutine npem

  subroutine npem_indrep(x,k,result,control,bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,.true.,ngrid,.false.)
  end subroutine npem_indrep

  subroutine npem_indrepbw(x,k,result,control,bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,.false.,ngrid,.false.)
  end subroutine npem_indrepbw

  subroutine npmsl(x,k,result,control,bandwidth,same_bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);logical,intent(in),optional::same_bandwidth;integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,same_bandwidth,ngrid,.false.)
  end subroutine npmsl

  subroutine spem(x,k,result,control,bandwidth,constant_bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);logical,intent(in),optional::constant_bandwidth;integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,constant_bandwidth,ngrid,.false.)
  end subroutine spem

  subroutine mvnpem(x,k,result,control,bandwidth,same_bandwidth,ngrid)
    real(dp),intent(in)::x(:,:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth(:);logical,intent(in),optional::same_bandwidth;integer,intent(in),optional::ngrid
    call product_kde_em(x,k,result,control,bandwidth,same_bandwidth,ngrid,.false.)
  end subroutine mvnpem

  subroutine spem_symloc(x,k,result,control,bandwidth,locations,ngrid)
    real(dp),intent(in)::x(:);integer,intent(in)::k
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::bandwidth,locations(:);integer,intent(in),optional::ngrid
    type(em_control)::ctl
    integer::n,m,i,j,iter,status
    real(dp),allocatable::post(:,:),lambda(:),mu(:),bw(:),history(:),lw(:),resid(:),dens(:),grid(:),densgrid(:,:),weights(:)
    real(dp)::ll,newll,ln,diff,h
    ctl=em_control();if(present(control))ctl=control;n=size(x);m=200;if(present(ngrid))m=ngrid
    allocate(post(n,k),lambda(k),mu(k),bw(k),history(ctl%max_iterations+1),lw(k),resid(n),dens(n),weights(n))
    call initialize_post(reshape(x,[n,1]),k,post);lambda=sum(post,dim=1)/real(n,dp)
    do j=1,k;mu(j)=sum(post(:,j)*x)/sum(post(:,j));end do
    if(present(locations))then;if(size(locations)==k)mu=locations;end if
    h=weighted_bandwidth(x,spread(1.0_dp,1,n));if(present(bandwidth))h=bandwidth;bw=h
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      lambda=sum(post,dim=1)/real(n,dp)
      do j=1,k;mu(j)=sum(post(:,j)*x)/max(sum(post(:,j)),tiny(1.0_dp));end do
      call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    call make_grid(x,m,grid);allocate(densgrid(k,m))
    do j=1,k;resid=x-mu(j);call wkde(resid,post(:,j),grid-mu(j),bw(j),densgrid(j,:),.true.);end do
    result%lambda=lambda;result%posterior=post;result%grid=grid;result%density=densgrid;result%bandwidth=bw;result%location=mu
    result%loglik=ll;result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=0
      do j=1,k
        resid=abs(x-mu(j));weights=post(:,j)
        call wkde(resid,weights,resid,bw(j),dens,.false.)
        post(:,j)=max(dens,tiny(1.0_dp))
      end do
      do i=1,n;do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+log(post(i,j));end do
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln;end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine spem_symloc

  subroutine spem_symloc_n01(x,result,control,location,bandwidth,ngrid)
    real(dp),intent(in)::x(:)
    type(semiparametric_result),intent(out)::result;type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::location,bandwidth;integer,intent(in),optional::ngrid
    type(em_control)::ctl
    integer::n,m,i,iter,status
    real(dp)::mu,h,ll,newll,ln,diff
    real(dp),allocatable::post(:,:),lambda(:),lw(:),res(:),dens(:),history(:),grid(:),dg(:,:),weights(:)
    ctl=em_control();if(present(control))ctl=control;n=size(x);m=200;if(present(ngrid))m=ngrid
    mu=sum(x)/real(n,dp)
    if(present(location))mu=location
    h=weighted_bandwidth(x,spread(1.0_dp,1,n))
    if(present(bandwidth))h=bandwidth
    allocate(post(n,2),lambda(2),lw(2),res(n),dens(n),history(ctl%max_iterations+1),weights(n));post=0.5_dp;lambda=0.5_dp
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      lambda=sum(post,dim=1)/real(n,dp);mu=sum(post(:,2)*x)/max(sum(post(:,2)),tiny(1.0_dp))
      call estep(newll);history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    call make_grid(x,m,grid);allocate(dg(2,m));dg(1,:)=normal_pdf(grid,0.0_dp,1.0_dp)
    res=x-mu;weights=post(:,2);call wkde(res,weights,grid-mu,h,dg(2,:),.true.)
    result%lambda=lambda;result%posterior=post;result%grid=grid;result%density=dg
    result%bandwidth=[1.0_dp,h];result%location=[0.0_dp,mu];result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0;result%status=merge(0,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      res=abs(x-mu);weights=post(:,2);call wkde(res,weights,res,h,dens,.false.);ll=0.0_dp;status=0
      do i=1,n;lw(1)=log(max(lambda(1),tiny(1.0_dp)))+log(max(normal_pdf(x(i),0.0_dp,1.0_dp),tiny(1.0_dp)))
        lw(2)=log(max(lambda(2),tiny(1.0_dp)))+log(max(dens(i),tiny(1.0_dp)))
        call normalize_logweights(lw,post(i,:),ln);ll=ll+ln;end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine spem_symloc_n01

  subroutine spregmix(y,x,k,regression_result,density_result,control,bandwidth,symmetric)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k
    type(regression_mixture_result),intent(out)::regression_result
    type(semiparametric_result),intent(out)::density_result
    type(em_control),intent(in),optional::control;real(dp),intent(in),optional::bandwidth
    logical,intent(in),optional::symmetric
    real(dp),allocatable::residual(:,:),bw(:),xa(:,:)
    integer::j,n
    call regmix_em(y,x,k,regression_result,control,.true.)
    n=size(y);allocate(residual(n,k),bw(k),xa(n,size(x,2)+1))
    xa(:,1)=1.0_dp;xa(:,2:)=x
    do j=1,k
      residual(:,j)=y-matmul(xa,regression_result%beta(:,j))
      bw(j)=weighted_bandwidth(residual(:,j),regression_result%posterior(:,j))
    end do
    if(present(bandwidth))bw=bandwidth
    call product_kde_em(residual,k,density_result,control,bw,.false.,200,symmetric)
  end subroutine spregmix
end module mixtools_semiparametric
