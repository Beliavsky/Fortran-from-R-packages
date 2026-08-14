module kriginv_criteria
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use kriginv_kinds, only : dp
  use kriginv_math, only : normal_pdf, normal_cdf, normal_quantile, bvn_cdf
  use kriginv_linalg, only : invert_spd
  use kriginv_model, only : krig_model, krig_prediction, predict_nobias_km, posterior_covariance
  implicit none
  private
  public :: excursion_probability, vorob_threshold
  public :: ranjan_optim, bichon_optim, tmse_optim, tsee_optim
  public :: predict_update_km_parallel, sur_optim_parallel, jn_optim_parallel
  public :: timse_optim_parallel, vorob_optim_parallel, vorobvol_optim_parallel
  public :: sur_optim_parallel2, jn_optim_parallel2, timse_optim_parallel2, vorob_optim_parallel2, vorobvol_optim_parallel2
  public :: compute_real_volume_constant
contains
  function excursion_probability(mean,sd,thresholds) result(p)
    real(dp), intent(in) :: mean(:),sd(:),thresholds(:)
    real(dp), allocatable :: p(:)
    integer :: i,j
    real(dp) :: z
    allocate(p(size(mean))); p=0.0_dp
    if(size(thresholds)==1) then
      do i=1,size(mean)
        if(sd(i)>0.0_dp) then
          p(i)=normal_cdf((mean(i)-thresholds(1))/sd(i))
        else
          p(i)=merge(1.0_dp,0.0_dp,mean(i)>thresholds(1))
        end if
      end do
    else
      do j=1,size(thresholds)
        do i=1,size(mean)
          if(sd(i)>0.0_dp) then
            z=(mean(i)-thresholds(j))/sd(i)
            p(i)=p(i)+merge(1.0_dp,-1.0_dp,mod(j,2)==1)*normal_cdf(z)
          else
            p(i)=p(i)+merge(1.0_dp,-1.0_dp,mod(j,2)==1)* &
                 merge(1.0_dp,0.0_dp,mean(i)>thresholds(j))
          end if
        end do
      end do
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function excursion_probability

  real(dp) function vorob_threshold(pn) result(alpha)
    real(dp), intent(in), optional :: pn(:)
    real(dp), allocatable :: s(:)
    real(dp) :: x,w
    integer :: n,k
    if(.not.present(pn)) then; alpha=0.5_dp; return; end if
    n=size(pn); if(n==0) then; alpha=0.5_dp; return; end if
    s=pn; call sort_real(s)
    x=(1.0_dp-sum(pn)/real(n,dp))*real(n,dp)
    ! R's original expression can index element zero.  Clamp the linear
    ! interpolation to the first/last order statistic instead.
    if(x<=1.0_dp) then
      alpha=s(1)
    else if(x>=real(n,dp)) then
      alpha=s(n)
    else
      k=floor(x); w=x-real(k,dp)
      alpha=(1.0_dp-w)*s(k)+w*s(k+1)
    end if
  end function vorob_threshold

  function ranjan_optim(x,model,threshold,alpha) result(crit)
    real(dp), intent(in) :: x(:,:),threshold
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: alpha
    real(dp), allocatable :: crit(:)
    type(krig_prediction) :: pr
    real(dp) :: al,t,tp,tm,g
    integer :: i
    al=1.0_dp; if(present(alpha)) al=alpha
    pr=predict_nobias_km(model,x,'UK',.false.); allocate(crit(size(x,1)))
    do i=1,size(x,1)
      if(pr%sd(i)<=0.0_dp) then; crit(i)=0.0_dp; cycle; end if
      t=(pr%mean(i)-threshold)/pr%sd(i); tp=t+al; tm=t-al
      g=(al*al-1.0_dp-t*t)*(normal_cdf(tp)-normal_cdf(tm))- &
        2.0_dp*t*(normal_pdf(tp)-normal_pdf(tm))+tp*normal_pdf(tp)-tm*normal_pdf(tm)
      if(ieee_is_nan(g)) g=0.0_dp
      crit(i)=g*pr%sd(i)**2
    end do
  end function ranjan_optim

  function bichon_optim(x,model,threshold,alpha) result(crit)
    real(dp), intent(in) :: x(:,:),threshold
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: alpha
    real(dp), allocatable :: crit(:)
    type(krig_prediction) :: pr
    real(dp) :: al,t,tp,tm,g
    integer :: i
    al=1.0_dp; if(present(alpha)) al=alpha
    pr=predict_nobias_km(model,x,'UK',.false.); allocate(crit(size(x,1)))
    do i=1,size(x,1)
      if(pr%sd(i)<=0.0_dp) then; crit(i)=0.0_dp; cycle; end if
      t=(pr%mean(i)-threshold)/pr%sd(i); tp=t+al; tm=t-al
      g=al*(normal_cdf(tp)-normal_cdf(tm))-t*(2.0_dp*normal_cdf(t)-normal_cdf(tp)-normal_cdf(tm))- &
        (2.0_dp*normal_pdf(t)-normal_pdf(tp)-normal_pdf(tm))
      crit(i)=g*pr%sd(i)
    end do
  end function bichon_optim

  function tmse_optim(x,model,thresholds,epsilon) result(crit)
    real(dp), intent(in) :: x(:,:),thresholds(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: epsilon
    real(dp), allocatable :: crit(:)
    type(krig_prediction) :: pr
    real(dp) :: eps,v,w
    integer :: i,j
    eps=0.0_dp; if(present(epsilon)) eps=epsilon
    pr=predict_nobias_km(model,x,'UK',.false.); allocate(crit(size(x,1))); crit=0.0_dp
    do i=1,size(x,1)
      v=pr%sd(i)**2+eps**2
      if(v<=0.0_dp) cycle
      w=0.0_dp
      do j=1,size(thresholds)
        w=w+normal_pdf((pr%mean(i)-thresholds(j))/sqrt(v))/sqrt(v)
      end do
      crit(i)=w*pr%sd(i)**2
    end do
  end function tmse_optim

  function tsee_optim(x,model,threshold) result(crit)
    real(dp), intent(in) :: x(:,:),threshold
    type(krig_model), intent(in) :: model
    real(dp), allocatable :: crit(:)
    type(krig_prediction) :: pr
    real(dp) :: t,a,b,sphi
    integer :: i
    pr=predict_nobias_km(model,x,'UK',.false.); allocate(crit(size(x,1))); crit=0.0_dp
    do i=1,size(x,1)
      if(pr%sd(i)<=0.0_dp) cycle
      t=(threshold-pr%mean(i))/pr%sd(i); sphi=pr%sd(i)*normal_pdf(t)
      a=(threshold-pr%mean(i))*normal_cdf(t)+sphi
      b=(pr%mean(i)-threshold)*normal_cdf(-t)+sphi
      crit(i)=a*b
    end do
  end function tsee_optim

  subroutine predict_update_km_parallel(newmean,sigma_r,newvalue,oldmean,oldsd,kn,newnoise,mean_out,sd_out,lambda,ok)
    real(dp), intent(in) :: newmean(:),sigma_r(:,:),newvalue(:),oldmean(:),oldsd(:),kn(:,:)
    real(dp), intent(in), optional :: newnoise(:)
    real(dp), allocatable, intent(out) :: mean_out(:),sd_out(:),lambda(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: s(:,:),sinv(:,:),diff(:),red(:)
    integer :: i
    s=sigma_r
    if(present(newnoise)) then
      if(size(newnoise)/=size(s,1)) then; ok=.false.; return; end if
      do i=1,size(s,1); s(i,i)=s(i,i)+max(0.0_dp,newnoise(i)); end do
    end if
    call invert_spd(s,sinv,ok)
    if(.not.ok) return
    lambda=matmul(kn,sinv); diff=newvalue-newmean
    mean_out=oldmean+matmul(lambda,diff)
    allocate(red(size(oldmean))); red=sum(lambda*kn,dim=2)
    sd_out=sqrt(max(0.0_dp,oldsd*oldsd-red))
  end subroutine predict_update_km_parallel

  real(dp) function sur_optim_parallel(xnew,integration_points,integration_weights,oldmean,oldsd,model,thresholds, &
                                       newnoise,current_sur) result(crit)
    real(dp), intent(in) :: xnew(:,:),integration_points(:,:),oldmean(:),oldsd(:),thresholds(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_sur
    type(krig_prediction) :: px
    real(dp), allocatable :: kn(:,:),nm(:),ns(:),lam(:,:),tmp(:)
    real(dp) :: rho,a
    integer :: i,j,k
    logical :: ok
    if(minimum_distance(model%x,xnew)<1.0e-5_dp .and. .not.present(newnoise)) then
      crit=1.01_dp*current_sur; return
    end if
    px=predict_nobias_km(model,xnew,'UK',.true.)
    kn=posterior_covariance(model,integration_points,xnew)
    call predict_update_km_parallel(px%mean,px%covariance,px%mean,oldmean,oldsd,kn,newnoise,nm,ns,lam,ok)
    if(.not.ok) then; crit=current_sur; return; end if
    allocate(tmp(size(oldmean))); tmp=0.0_dp
    if(size(thresholds)==1) then
      do i=1,size(oldmean)
        if(oldsd(i)<=0.0_dp) cycle
        a=(oldmean(i)-thresholds(1))/oldsd(i)
        rho=-(oldsd(i)**2-ns(i)**2)/max(oldsd(i)**2,tiny(1.0_dp))
        rho=max(-1.0_dp,min(1.0_dp,rho))
        tmp(i)=bvn_cdf(a,-a,rho)
      end do
    else
      do i=1,size(oldmean)
        if(oldsd(i)<=0.0_dp) cycle
        rho=(oldsd(i)**2-ns(i)**2)/max(oldsd(i)**2,tiny(1.0_dp))
        rho=max(-1.0_dp,min(1.0_dp,rho))
        do j=1,size(thresholds)
          a=(oldmean(i)-thresholds(j))/oldsd(i)
          tmp(i)=tmp(i)-parity_sign(j)*normal_cdf(a)
          do k=1,size(thresholds)
            tmp(i)=tmp(i)-parity_sign(j+k)*bvn_cdf(a,(oldmean(i)-thresholds(k))/oldsd(i),rho)
          end do
        end do
      end do
    end if
    crit=weighted_average(tmp,integration_weights)
  end function sur_optim_parallel

  real(dp) function timse_optim_parallel(xnew,integration_points,integration_weights,oldmean,oldsd,model, &
                                         weight,newnoise,current_timse) result(crit)
    real(dp), intent(in) :: xnew(:,:),integration_points(:,:),oldmean(:),oldsd(:),weight(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_timse
    type(krig_prediction) :: px
    real(dp), allocatable :: kn(:,:),nm(:),ns(:),lam(:,:),v(:)
    logical :: ok
    if(minimum_distance(model%x,xnew)<1.0e-5_dp .and. .not.present(newnoise)) then
      crit=1.01_dp*current_timse; return
    end if
    px=predict_nobias_km(model,xnew,'UK',.true.); kn=posterior_covariance(model,integration_points,xnew)
    call predict_update_km_parallel(px%mean,px%covariance,px%mean,oldmean,oldsd,kn,newnoise,nm,ns,lam,ok)
    if(.not.ok) then; crit=current_timse; return; end if
    v=weight*ns*ns; crit=weighted_average(v,integration_weights)
  end function timse_optim_parallel

  real(dp) function vorob_optim_parallel(xnew,integration_points,integration_weights,oldmean,oldsd,model,threshold, &
                                         alpha,penalisation,type_ex,newnoise,current_vorob) result(crit)
    real(dp), intent(in) :: xnew(:,:),integration_points(:,:),oldmean(:),oldsd(:),threshold,alpha
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:),penalisation
    character(len=*), intent(in), optional :: type_ex
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_vorob
    type(krig_prediction) :: px
    real(dp), allocatable :: kn(:,:),nm(:),ns(:),lam(:,:),res(:)
    real(dp) :: pen,a,c,arg1,arg2,arg3,q
    integer :: i
    logical :: ok,above
    pen=1.0_dp; if(present(penalisation)) pen=penalisation
    above=.true.; if(present(type_ex)) above=(type_ex(1:1)=='>')
    if(minimum_distance(model%x,xnew)<1.0e-5_dp .and. .not.present(newnoise)) then
      crit=current_vorob+0.01_dp; return
    end if
    px=predict_nobias_km(model,xnew,'UK',.true.); kn=posterior_covariance(model,integration_points,xnew)
    call predict_update_km_parallel(px%mean,px%covariance,px%mean,oldmean,oldsd,kn,newnoise,nm,ns,lam,ok)
    if(.not.ok) then; crit=current_vorob; return; end if
    allocate(res(size(oldmean))); q=normal_quantile(alpha)
    do i=1,size(oldmean)
      if(oldsd(i)<=0.0_dp .or. ns(i)<=0.0_dp) then; res(i)=0.0_dp; cycle; end if
      a=(oldmean(i)-threshold)/ns(i); c=oldsd(i)**2/ns(i)**2
      if(c<=1.0_dp+1.0e-14_dp) then; res(i)=current_vorob; cycle; end if
      arg1=(oldmean(i)-threshold)/oldsd(i); arg3=-sqrt(max(0.0_dp,1.0_dp-1.0_dp/c))
      if(above) then
        arg2=(q-a)/sqrt(c-1.0_dp)
        res(i)=bvn_cdf(arg1,arg2,arg3)-pen*bvn_cdf(arg1,-arg2,-arg3)+pen*normal_cdf(-arg2)
      else
        arg2=(q+a)/sqrt(c-1.0_dp)
        res(i)=bvn_cdf(-arg1,arg2,arg3)-pen*bvn_cdf(-arg1,-arg2,-arg3)+ &
               pen*normal_cdf((-q-a)/sqrt(c-1.0_dp))
      end if
    end do
    crit=weighted_average(res,integration_weights)
  end function vorob_optim_parallel

  real(dp) function vorobvol_optim_parallel(xnew,integration_points,integration_weights,oldmean,oldsd,model,threshold, &
                                            alpha,type_ex,newnoise,current_crit) result(crit)
    real(dp), intent(in) :: xnew(:,:),integration_points(:,:),oldmean(:),oldsd(:),threshold,alpha
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    character(len=*), intent(in), optional :: type_ex
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_crit
    type(krig_prediction) :: px
    real(dp), allocatable :: kn(:,:),nm(:),ns(:),lam(:,:),res(:)
    real(dp) :: a,c,arg,q
    integer :: i
    logical :: ok,above
    above=.true.; if(present(type_ex)) above=(type_ex(1:1)=='>')
    if(minimum_distance(model%x,xnew)<1.0e-5_dp .and. .not.present(newnoise)) then
      crit=current_crit+0.01_dp; return
    end if
    px=predict_nobias_km(model,xnew,'UK',.true.); kn=posterior_covariance(model,integration_points,xnew)
    call predict_update_km_parallel(px%mean,px%covariance,px%mean,oldmean,oldsd,kn,newnoise,nm,ns,lam,ok)
    if(.not.ok) then; crit=current_crit; return; end if
    allocate(res(size(oldmean))); q=normal_quantile(alpha)
    do i=1,size(oldmean)
      if(ns(i)<=0.0_dp .or. oldsd(i)<=0.0_dp) then; res(i)=0.0_dp; cycle; end if
      a=(oldmean(i)-threshold)/ns(i); c=oldsd(i)**2/ns(i)**2
      if(c<=1.0_dp+1.0e-14_dp) then; res(i)=0.0_dp; cycle; end if
      if(above) then; arg=(q-a)/sqrt(c-1.0_dp)
      else; arg=(q+a)/sqrt(c-1.0_dp); end if
      res(i)=normal_cdf(-arg)
    end do
    crit=weighted_average(res,integration_weights)
  end function vorobvol_optim_parallel

  real(dp) function jn_optim_parallel(xnew,integration_points,integration_weights,oldmean,oldsd,model,thresholds, &
                                      newnoise,current_sur) result(crit)
    real(dp), intent(in) :: xnew(:,:),integration_points(:,:),oldmean(:),oldsd(:),thresholds(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_sur
    type(krig_prediction) :: px
    real(dp), allocatable :: kn(:,:),nm(:),ns(:),lam(:,:),tmp(:),sigma_l(:,:)
    real(dp) :: rho,a,b,w
    integer :: i,j,u,v,m,np
    logical :: ok,paired
    if(minimum_distance(model%x,xnew)<1.0e-5_dp .and. .not.present(newnoise)) then
      crit=current_sur; return
    end if
    px=predict_nobias_km(model,xnew,'UK',.true.); kn=posterior_covariance(model,integration_points,xnew)
    call predict_update_km_parallel(px%mean,px%covariance,px%mean,oldmean,oldsd,kn,newnoise,nm,ns,lam,ok)
    if(.not.ok) then; crit=current_sur; return; end if
    sigma_l=matmul(matmul(lam,px%covariance),transpose(lam))
    paired=present(integration_weights)
    if(paired) paired=(size(integration_points,1)==2*size(integration_weights))
    if(paired) then
      m=size(integration_weights); allocate(tmp(m)); tmp=0.0_dp
      do i=1,m
        if(oldsd(i)<=0.0_dp .or. oldsd(m+i)<=0.0_dp) cycle
        rho=sigma_l(i,m+i)/(oldsd(i)*oldsd(m+i)); rho=max(-1.0_dp,min(1.0_dp,rho))
        do u=1,size(thresholds)
          a=(oldmean(i)-thresholds(u))/oldsd(i)
          do v=1,size(thresholds)
            b=(oldmean(m+i)-thresholds(v))/oldsd(m+i)
            tmp(i)=tmp(i)+parity_sign(u+v)*bvn_cdf(a,b,rho)
          end do
        end do
      end do
      crit=-sum(tmp*integration_weights)
    else
      np=size(integration_points,1); crit=0.0_dp
      do i=1,np
        if(oldsd(i)<=0.0_dp) cycle
        do j=1,np
          if(oldsd(j)<=0.0_dp) cycle
          rho=sigma_l(i,j)/(oldsd(i)*oldsd(j)); rho=max(-1.0_dp,min(1.0_dp,rho)); w=1.0_dp
          if(present(integration_weights)) w=integration_weights(i)*integration_weights(j)
          do u=1,size(thresholds)
            a=(oldmean(i)-thresholds(u))/oldsd(i)
            do v=1,size(thresholds)
              b=(oldmean(j)-thresholds(v))/oldsd(j)
              crit=crit-w*parity_sign(u+v)*bvn_cdf(a,b,rho)
            end do
          end do
        end do
      end do
      if(.not.present(integration_weights)) crit=crit/real(np*np,dp)
    end if
  end function jn_optim_parallel


  real(dp) function sur_optim_parallel2(x,other_points,integration_points,integration_weights,oldmean,oldsd,model, &
                                        thresholds,newnoise,current_sur) result(crit)
    real(dp), intent(in) :: x(:),other_points(:,:),integration_points(:,:),oldmean(:),oldsd(:),thresholds(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_sur
    real(dp), allocatable :: allx(:,:)
    allx=append_candidate(x,other_points)
    if(present(integration_weights)) then
      crit=sur_optim_parallel(allx,integration_points,integration_weights,oldmean,oldsd,model,thresholds,newnoise,current_sur)
    else
      crit=sur_optim_parallel(allx,integration_points,oldmean=oldmean,oldsd=oldsd,model=model,thresholds=thresholds, &
                              newnoise=newnoise,current_sur=current_sur)
    end if
  end function sur_optim_parallel2

  real(dp) function jn_optim_parallel2(x,other_points,integration_points,integration_weights,oldmean,oldsd,model, &
                                       thresholds,newnoise,current_sur) result(crit)
    real(dp), intent(in) :: x(:),other_points(:,:),integration_points(:,:),oldmean(:),oldsd(:),thresholds(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_sur
    real(dp), allocatable :: allx(:,:)
    allx=append_candidate(x,other_points)
    if(present(integration_weights)) then
      crit=jn_optim_parallel(allx,integration_points,integration_weights,oldmean,oldsd,model,thresholds,newnoise,current_sur)
    else
      crit=jn_optim_parallel(allx,integration_points,oldmean=oldmean,oldsd=oldsd,model=model,thresholds=thresholds, &
                             newnoise=newnoise,current_sur=current_sur)
    end if
  end function jn_optim_parallel2

  real(dp) function timse_optim_parallel2(x,other_points,integration_points,integration_weights,oldmean,oldsd,model, &
                                          weight,newnoise,current_timse) result(crit)
    real(dp), intent(in) :: x(:),other_points(:,:),integration_points(:,:),oldmean(:),oldsd(:),weight(:)
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_timse
    real(dp), allocatable :: allx(:,:)
    allx=append_candidate(x,other_points)
    if(present(integration_weights)) then
      crit=timse_optim_parallel(allx,integration_points,integration_weights,oldmean,oldsd,model,weight,newnoise,current_timse)
    else
      crit=timse_optim_parallel(allx,integration_points,oldmean=oldmean,oldsd=oldsd,model=model,weight=weight, &
                                newnoise=newnoise,current_timse=current_timse)
    end if
  end function timse_optim_parallel2

  real(dp) function vorob_optim_parallel2(x,other_points,integration_points,integration_weights,oldmean,oldsd,model, &
                                          threshold,alpha,penalisation,type_ex,newnoise,current_vorob) result(crit)
    real(dp), intent(in) :: x(:),other_points(:,:),integration_points(:,:),oldmean(:),oldsd(:),threshold,alpha
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:),penalisation
    character(len=*), intent(in), optional :: type_ex
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_vorob
    real(dp), allocatable :: allx(:,:)
    allx=append_candidate(x,other_points)
    if(present(integration_weights)) then
      crit=vorob_optim_parallel(allx,integration_points,integration_weights,oldmean,oldsd,model,threshold,alpha, &
                                penalisation,type_ex,newnoise,current_vorob)
    else
      crit=vorob_optim_parallel(allx,integration_points,oldmean=oldmean,oldsd=oldsd,model=model,threshold=threshold, &
                                alpha=alpha,penalisation=penalisation,type_ex=type_ex,newnoise=newnoise,current_vorob=current_vorob)
    end if
  end function vorob_optim_parallel2

  real(dp) function vorobvol_optim_parallel2(x,other_points,integration_points,integration_weights,oldmean,oldsd,model, &
                                             threshold,alpha,type_ex,newnoise,current_crit) result(crit)
    real(dp), intent(in) :: x(:),other_points(:,:),integration_points(:,:),oldmean(:),oldsd(:),threshold,alpha
    real(dp), intent(in), optional :: integration_weights(:),newnoise(:)
    character(len=*), intent(in), optional :: type_ex
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: current_crit
    real(dp), allocatable :: allx(:,:)
    allx=append_candidate(x,other_points)
    if(present(integration_weights)) then
      crit=vorobvol_optim_parallel(allx,integration_points,integration_weights,oldmean,oldsd,model,threshold,alpha, &
                                   type_ex,newnoise,current_crit)
    else
      crit=vorobvol_optim_parallel(allx,integration_points,oldmean=oldmean,oldsd=oldsd,model=model,threshold=threshold, &
                                   alpha=alpha,type_ex=type_ex,newnoise=newnoise,current_crit=current_crit)
    end if
  end function vorobvol_optim_parallel2

  function append_candidate(x,other) result(allx)
    real(dp), intent(in) :: x(:),other(:,:)
    real(dp), allocatable :: allx(:,:)
    integer :: n
    n=size(other,1); allocate(allx(n+1,size(x)))
    if(n>0) allx(1:n,:)=other
    allx(n+1,:)=x
  end function append_candidate

  real(dp) function compute_real_volume_constant(model,points,threshold,weights) result(res)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: points(:,:),threshold
    real(dp), intent(in), optional :: weights(:)
    type(krig_prediction) :: pr
    real(dp) :: a,b,rho,w
    integer :: i,j,n
    pr=predict_nobias_km(model,points,'UK',.true.); n=size(points,1); res=0.0_dp
    do i=1,n
      if(pr%sd(i)<=0.0_dp) cycle
      a=(threshold-pr%mean(i))/pr%sd(i)
      do j=1,n
        if(pr%sd(j)<=0.0_dp) cycle
        b=(threshold-pr%mean(j))/pr%sd(j)
        rho=pr%covariance(i,j)/(pr%sd(i)*pr%sd(j)); rho=max(-1.0_dp,min(1.0_dp,rho)); w=1.0_dp
        if(present(weights)) w=weights(i)*weights(j)
        res=res+w*(1.0_dp-normal_cdf(a)-normal_cdf(b)+bvn_cdf(a,b,rho))
      end do
    end do
    if(.not.present(weights)) res=res/real(n*n,dp)
  end function compute_real_volume_constant

  real(dp) function weighted_average(x,w) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: w(:)
    if(present(w)) then; v=sum(x*w)
    else; v=sum(x)/real(max(1,size(x)),dp); end if
  end function weighted_average

  real(dp) function minimum_distance(oldx,newx) result(v)
    real(dp), intent(in) :: oldx(:,:),newx(:,:)
    integer :: i,j
    v=huge(1.0_dp)
    do i=1,size(newx,1)
      do j=1,size(oldx,1); v=min(v,sqrt(sum((newx(i,:)-oldx(j,:))**2))); end do
      do j=1,i-1; v=min(v,sqrt(sum((newx(i,:)-newx(j,:))**2))); end do
    end do
  end function minimum_distance

  real(dp) function parity_sign(i) result(s)
    integer, intent(in) :: i
    if(mod(i,2)==0) then; s=1.0_dp; else; s=-1.0_dp; end if
  end function parity_sign

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: t
    do i=2,size(x)
      t=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=t) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=t
    end do
  end subroutine sort_real
end module kriginv_criteria
