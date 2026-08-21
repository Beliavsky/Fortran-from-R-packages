module relsurv_years
  use relsurv_kinds, only : dp
  implicit none
  private

  type, public :: years_result
    real(dp), allocatable :: time(:), estimate(:), standard_error(:), lower(:), upper(:)
    real(dp), allocatable :: observed_prob(:), population_prob(:)
    real(dp), allocatable :: observed_prob_se(:), population_prob_se(:)
    real(dp), allocatable :: observed_area(:), population_area(:)
    real(dp), allocatable :: observed_area_se(:), population_area_se(:)
  end type years_result

  type, public :: yl2013_result
    real(dp), allocatable :: time(:), excess_failure(:), excess_area(:)
    real(dp), allocatable :: estimate(:), standard_error(:), lower(:), upper(:)
    real(dp), allocatable :: excess_failure_se(:), excess_area_se(:)
  end type yl2013_result

  public :: years_difference, years_yl2013, years_yl2017
  public :: greenwood_area_variance, bootstrap_column_variance

contains

  subroutine years_difference(time,survival,n_risk,n_event,pop_failure,result,scale,conf_int, &
      boot_failure,boot_pop_failure)
    real(dp),intent(in)::time(:),survival(:),n_risk(:),n_event(:),pop_failure(:)
    type(years_result),intent(out)::result
    real(dp),intent(in),optional::scale,conf_int
    real(dp),intent(in),optional::boot_failure(:,:),boot_pop_failure(:,:)
    integer::n
    real(dp)::sc,ci,z
    real(dp),allocatable::varea(:),vprob(:),barea(:,:),bparea(:,:),tmp(:)
    n=size(time)
    if(any([size(survival),size(n_risk),size(n_event),size(pop_failure)]/=n))error stop 'years_difference: shape'
    sc=365.241_dp;if(present(scale))sc=scale
    ci=0.95_dp;if(present(conf_int))ci=conf_int
    z=normal_quantile(0.5_dp+ci/2.0_dp)
    allocate(result%time(n),result%estimate(n),result%standard_error(n),result%lower(n),result%upper(n))
    allocate(result%observed_prob(n),result%population_prob(n),result%observed_prob_se(n),result%population_prob_se(n))
    allocate(result%observed_area(n),result%population_area(n),result%observed_area_se(n),result%population_area_se(n))
    result%time=time;result%observed_prob=1.0_dp-survival;result%population_prob=pop_failure
    call left_area(time,result%observed_prob,sc,result%observed_area)
    call left_area(time,result%population_prob,sc,result%population_area)
    result%estimate=result%observed_area-result%population_area
    call greenwood_prob_se(survival,n_risk,n_event,result%observed_prob_se)
    call greenwood_area_variance(time,survival,n_risk,n_event,varea,sc)
    result%observed_area_se=sqrt(max(varea,0.0_dp));result%population_prob_se=0.0_dp;result%population_area_se=0.0_dp
    result%standard_error=result%observed_area_se

    if(present(boot_failure).neqv.present(boot_pop_failure))error stop 'years_difference: provide both bootstrap curve arrays'
    if(present(boot_failure))then
      if(size(boot_failure,2)/=n.or.size(boot_pop_failure,2)/=n.or.size(boot_failure,1)/=size(boot_pop_failure,1)) &
        error stop 'years_difference: bootstrap shape'
      call bootstrap_column_variance(boot_failure,vprob);result%observed_prob_se=sqrt(max(vprob,0.0_dp))
      call bootstrap_column_variance(boot_pop_failure,vprob);result%population_prob_se=sqrt(max(vprob,0.0_dp))
      allocate(barea(size(boot_failure,1),n),bparea(size(boot_failure,1),n))
      call bootstrap_areas(time,boot_failure,sc,barea);call bootstrap_areas(time,boot_pop_failure,sc,bparea)
      call bootstrap_column_variance(barea,varea);result%observed_area_se=sqrt(max(varea,0.0_dp))
      call bootstrap_column_variance(bparea,varea);result%population_area_se=sqrt(max(varea,0.0_dp))
      call bootstrap_difference_variance(barea,bparea,varea);result%standard_error=sqrt(max(varea,0.0_dp))
    end if
    result%lower=result%estimate-z*result%standard_error;result%upper=result%estimate+z*result%standard_error
  end subroutine years_difference

  subroutine years_yl2017(time,observed_failure,population_failure,result,scale,conf_int, &
      boot_observed,boot_population)
    real(dp),intent(in)::time(:),observed_failure(:),population_failure(:)
    type(years_result),intent(out)::result
    real(dp),intent(in),optional::scale,conf_int
    real(dp),intent(in),optional::boot_observed(:,:),boot_population(:,:)
    real(dp),allocatable::surv(:),nr(:),ne(:)
    allocate(surv(size(time)),nr(size(time)),ne(size(time)))
    surv=1.0_dp-observed_failure;nr=huge(1.0_dp);ne=0.0_dp
    call years_difference(time,surv,nr,ne,population_failure,result,scale,conf_int,boot_observed,boot_population)
  end subroutine years_yl2017

  subroutine years_yl2013(time,survival,n_risk,n_event,pop_hazard,result,scale,conf_int,boot_excess_failure)
    real(dp),intent(in)::time(:),survival(:),n_risk(:),n_event(:),pop_hazard(:)
    type(yl2013_result),intent(out)::result
    real(dp),intent(in),optional::scale,conf_int
    real(dp),intent(in),optional::boot_excess_failure(:,:)
    integer::n,i
    real(dp)::sc,ci,z
    real(dp),allocatable::v(:),ba(:,:)
    n=size(time);if(any([size(survival),size(n_risk),size(n_event),size(pop_hazard)]/=n))error stop 'years_yl2013: shape'
    sc=365.241_dp;if(present(scale))sc=scale;ci=0.95_dp;if(present(conf_int))ci=conf_int
    z=normal_quantile(0.5_dp+ci/2.0_dp)
    allocate(result%time(n),result%excess_failure(n),result%excess_area(n),result%estimate(n))
    allocate(result%standard_error(n),result%lower(n),result%upper(n),result%excess_failure_se(n),result%excess_area_se(n))
    result%time=time;result%excess_failure=0.0_dp
    do i=1,n
      if(i==1)then
        result%excess_failure(i)=survival(i)*(n_event(i)/max(n_risk(i),tiny(1.0_dp))-pop_hazard(i))
      else
        result%excess_failure(i)=result%excess_failure(i-1)+ &
          survival(i)*(n_event(i)/max(n_risk(i),tiny(1.0_dp))-pop_hazard(i))
      end if
    end do
    call left_area(time,result%excess_failure,sc,result%excess_area)
    ! Upstream YL2013 estimate uses cumsum(F_E * c(0,diff(time))), i.e. right value
    result%estimate=0.0_dp
    do i=2,n;result%estimate(i)=result%estimate(i-1)+(time(i)-time(i-1))*result%excess_failure(i)/sc;end do
    result%standard_error=0.0_dp;result%excess_failure_se=0.0_dp;result%excess_area_se=0.0_dp
    if(present(boot_excess_failure))then
      if(size(boot_excess_failure,2)/=n)error stop 'years_yl2013: bootstrap shape'
      call bootstrap_column_variance(boot_excess_failure,v);result%excess_failure_se=sqrt(max(v,0.0_dp))
      allocate(ba(size(boot_excess_failure,1),n));call bootstrap_areas(time,boot_excess_failure,sc,ba)
      call bootstrap_column_variance(ba,v);result%excess_area_se=sqrt(max(v,0.0_dp))
      call bootstrap_yl2013_estimates(time,boot_excess_failure,sc,ba)
      call bootstrap_column_variance(ba,v);result%standard_error=sqrt(max(v,0.0_dp))
    end if
    result%lower=result%estimate-z*result%standard_error;result%upper=result%estimate+z*result%standard_error
  end subroutine years_yl2013

  subroutine greenwood_prob_se(survival,nrisk,nevent,se)
    real(dp),intent(in)::survival(:),nrisk(:),nevent(:)
    real(dp),intent(out)::se(:)
    real(dp)::acc,den
    integer::i,lastgood
    acc=0.0_dp;lastgood=0
    do i=1,size(survival)
      den=nrisk(i)*(nrisk(i)-nevent(i))
      if(den>0.0_dp)acc=acc+nevent(i)/den
      se(i)=sqrt(max(0.0_dp,survival(i)*survival(i)*acc))
      if(survival(i)>0.0_dp)lastgood=i
    end do
    if(lastgood>0.and.lastgood<size(se))se(lastgood+1:)=se(lastgood)
  end subroutine greenwood_prob_se

  subroutine greenwood_area_variance(time,survival,nrisk,nevent,var,scale)
    real(dp),intent(in)::time(:),survival(:),nrisk(:),nevent(:)
    real(dp),allocatable,intent(out)::var(:)
    real(dp),intent(in),optional::scale
    real(dp)::sc,dt,den,s
    real(dp),allocatable::auc(:),leftsurv(:)
    integer::i,j,n
    n=size(time);sc=365.241_dp;if(present(scale))sc=scale
    allocate(var(n),auc(n),leftsurv(n));leftsurv(1)=1.0_dp;if(n>1)leftsurv(2:)=survival(:n-1)
    auc=0.0_dp
    if(n>0)auc(1)=0.0_dp
    do i=2,n
      dt=time(i)-time(i-1)
      auc(i)=auc(i-1)+dt*leftsurv(i)/sc
    end do
    do i=1,n
      s=0.0_dp
      do j=1,i
        den=nrisk(j)*(nrisk(j)-nevent(j))
        if(den>0.0_dp)s=s+(auc(i)-auc(j))**2*nevent(j)/den
      end do
      var(i)=s
    end do
    if(n>1)then
      do i=2,n
        if(var(i)/=var(i))var(i)=var(i-1)
      end do
    end if
  end subroutine greenwood_area_variance

  subroutine left_area(time,failure,scale,area)
    real(dp),intent(in)::time(:),failure(:),scale
    real(dp),intent(out)::area(:)
    real(dp)::dt,prev
    integer::i
    area=0.0_dp;prev=0.0_dp
    if(size(time)>0)area(1)=0.0_dp
    do i=2,size(time)
      dt=time(i)-time(i-1);area(i)=area(i-1)+dt*failure(i-1)/scale
    end do
  end subroutine left_area

  subroutine bootstrap_areas(time,curves,scale,areas)
    real(dp),intent(in)::time(:),curves(:,:),scale
    real(dp),intent(out)::areas(:,:)
    integer::b
    do b=1,size(curves,1);call left_area(time,curves(b,:),scale,areas(b,:));end do
  end subroutine bootstrap_areas

  subroutine bootstrap_yl2013_estimates(time,curves,scale,est)
    real(dp),intent(in)::time(:),curves(:,:),scale
    real(dp),intent(out)::est(:,:)
    integer::b,i
    est=0.0_dp
    do b=1,size(curves,1)
      do i=2,size(time);est(b,i)=est(b,i-1)+(time(i)-time(i-1))*curves(b,i)/scale;end do
    end do
  end subroutine bootstrap_yl2013_estimates

  subroutine bootstrap_column_variance(x,var)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::var(:)
    integer::j,n
    real(dp)::m
    allocate(var(size(x,2)));n=size(x,1)
    do j=1,size(x,2)
      if(n<=1)then;var(j)=0.0_dp
      else;m=sum(x(:,j))/real(n,dp);var(j)=sum((x(:,j)-m)**2)/real(n-1,dp);end if
    end do
  end subroutine bootstrap_column_variance

  subroutine bootstrap_difference_variance(a,b,var)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable,intent(out)::var(:)
    real(dp),allocatable::d(:,:)
    allocate(d(size(a,1),size(a,2)));d=a-b;call bootstrap_column_variance(d,var)
  end subroutine bootstrap_difference_variance

  function normal_quantile(p) result(x)
    real(dp),intent(in)::p
    real(dp)::x,q,r
    real(dp),parameter::a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp
    real(dp),parameter::a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp
    real(dp),parameter::a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
    real(dp),parameter::b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp
    real(dp),parameter::b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp
    real(dp),parameter::c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp
    real(dp),parameter::c3=-2.400758277161838_dp,c4=-2.549732539343734_dp
    real(dp),parameter::c5=4.374664141464968_dp,c6=2.938163982698783_dp
    real(dp),parameter::d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp
    real(dp),parameter::d3=2.445134137142996_dp,d4=3.754408661907416_dp
    if(p<=0.0_dp)then;x=-huge(1.0_dp);return;end if
    if(p>=1.0_dp)then;x=huge(1.0_dp);return;end if
    if(p<0.02425_dp)then
      q=sqrt(-2.0_dp*log(p));x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if(p>0.97575_dp)then
      q=sqrt(-2.0_dp*log(1.0_dp-p));x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else
      q=p-0.5_dp;r=q*q
      x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    end if
  end function normal_quantile

end module relsurv_years
