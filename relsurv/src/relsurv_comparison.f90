module relsurv_comparison
  use relsurv_kinds, only : dp
  use relsurv_ratetable, only : ratetable_type, net_summary_type, netfast_summary, netwei_summary, pystep
  use relsurv_linalg, only : solve_linear
  use relsurv_stats, only : chi_square_sf
  implicit none
  private
  type, public :: cmp_rel_result
    real(dp), allocatable :: time(:), disease(:), population(:), var_disease(:), var_population(:)
    real(dp), allocatable :: area_disease_increment(:), area_population_increment(:)
    real(dp) :: area_disease=0.0_dp, area_population=0.0_dp
  end type cmp_rel_result
  type, public :: rsdiff_result
    real(dp), allocatable :: statistic_vector(:), covariance(:,:)
    real(dp) :: chisq=0.0_dp, p_value=1.0_dp
    integer :: df=0
  end type rsdiff_result
  public :: cmp_rel, rsdiff, nessie_expected
contains
  subroutine cmp_rel(tab,x,y,status,times,result,ystart,scale)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),y(:),times(:)
    integer,intent(in)::status(:)
    type(cmp_rel_result),intent(out)::result
    real(dp),intent(in),optional::ystart(:),scale
    real(dp),allocatable::ys(:),si(:),sitt(:),dlp(:),dlo(:),dle(:),sigma(:),sigmap(:),so(:),soprev(:)
    real(dp)::data(tab%ndim),time0,thiscell,remain,dt,wt,hazard,sc
    integer::i,j,k,kt,idx,idx2,n
    type(net_summary_type)::s
    n=size(y); allocate(ys(n)); ys=0.0_dp; if(present(ystart))ys=ystart
    ! The upstream cmpfast population/event sums are identical to the exact net accumulator for these fields.
    call netwei_summary(tab,x,y,status,times,s)
    allocate(result%time(size(times)),result%disease(size(times)),result%population(size(times)), &
      result%var_disease(size(times)),result%var_population(size(times)), &
      result%area_disease_increment(size(times)),result%area_population_increment(size(times)))
    allocate(si(n),sitt(n),dlp(size(times)),dlo(size(times)),dle(size(times)),sigma(size(times)), &
      sigmap(size(times)),so(size(times)),soprev(size(times)))
    ! Recompute yidli respecting late entry exactly as cmpfast.
    s%yidli=0.0_dp; s%yi=0.0_dp; s%dni=0.0_dp; si=1.0_dp; time0=0.0_dp
    do j=1,size(times)
      thiscell=times(j)-time0
      do i=1,n
        if(y(i)>=times(j)) then
          data=x(i,:)
          do k=1,tab%ndim; if(tab%factor(k)/=1)data(k)=data(k)+time0; end do
          remain=thiscell; hazard=0.0_dp
          do while(remain>100.0_dp*epsilon(1.0_dp))
            call pystep(tab,data,remain,1,dt,idx,idx2,wt); if(dt<=0.0_dp)exit
            if(wt<1.0_dp)then
              hazard=hazard+dt*(wt*tab%rate(idx)+(1.0_dp-wt)*tab%rate(idx2))
            else
              hazard=hazard+dt*tab%rate(idx)
            end if
            do k=1,tab%ndim; if(tab%factor(k)/=1)data(k)=data(k)+dt; end do
            remain=remain-dt
          end do
          sitt(i)=si(i); si(i)=si(i)*exp(-hazard)
          if(ys(i)<times(j))then
            s%yidli(j)=s%yidli(j)+hazard; s%yi(j)=s%yi(j)+1.0_dp
            if(y(i)==times(j))s%dni(j)=s%dni(j)+real(status(i),dp)
          end if
        end if
      end do
      time0=times(j)
    end do
    result%disease=0.0_dp; result%population=0.0_dp; result%var_disease=0.0_dp; result%var_population=0.0_dp
    so=1.0_dp; soprev=1.0_dp
    do j=1,size(times)
      if(s%yi(j)>0.0_dp)then
        dlp(j)=s%yidli(j)/s%yi(j); dlo(j)=s%dni(j)/s%yi(j)
        sigma(j)=s%dni(j)/(s%yi(j)*s%yi(j)); sigmap(j)=s%yidli(j)/(s%yi(j)*s%yi(j))
      else
        dlp(j)=0.0_dp; dlo(j)=0.0_dp; sigma(j)=0.0_dp; sigmap(j)=0.0_dp
      end if
      dle(j)=dlo(j)-dlp(j)
      if(j==1)then; soprev(j)=1.0_dp; so(j)=1.0_dp-dlo(j)
      else; soprev(j)=so(j-1); so(j)=so(j-1)*(1.0_dp-dlo(j)); end if
      if(j==1)then
        result%disease(j)=soprev(j)*dle(j); result%population(j)=soprev(j)*dlp(j)
      else
        result%disease(j)=result%disease(j-1)+soprev(j)*dle(j)
        result%population(j)=result%population(j-1)+soprev(j)*dlp(j)
      end if
      do kt=1,j
        if(so(kt)>tiny(1.0_dp))then
          result%var_disease(j)=result%var_disease(j)+so(kt)**2* &
            (1.0_dp-(result%disease(j)-result%disease(kt))/so(kt))**2*sigma(kt)
        end if
        result%var_population(j)=result%var_population(j)+ &
          (result%population(j)-result%population(kt))**2*sigma(kt)
      end do
      result%area_disease_increment(j)=this_interval(times,j)*result%disease(j)
      result%area_population_increment(j)=this_interval(times,j)*result%population(j)
    end do
    result%time=times; sc=365.241_dp; if(present(scale))sc=scale
    result%area_disease=sum(result%area_disease_increment)/sc
    result%area_population=sum(result%area_population_increment)/sc
  contains
    pure function this_interval(t,j) result(dt0)
      real(dp),intent(in)::t(:); integer,intent(in)::j; real(dp)::dt0
      if(j==1)then; dt0=t(1); else; dt0=t(j)-t(j-1); end if
    end function
  end subroutine cmp_rel

  subroutine rsdiff(tab,x,y,status,group,times,result,stratum,precision)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),y(:),times(:)
    integer,intent(in)::status(:),group(:)
    type(rsdiff_result),intent(out)::result
    integer,intent(in),optional::stratum(:)
    real(dp),intent(in),optional::precision
    integer::gmax,smax,g,s,t,n,ng,ns,i,j,k
    integer,allocatable::str(:)
    real(dp),allocatable::wr(:,:,:),we(:,:,:),dn2(:,:,:),totalr(:,:),totale(:,:),z(:),cov(:,:),v(:),sol(:)
    real(dp),allocatable::xx(:,:),yy(:),ys(:)
    integer,allocatable::ev(:)
    type(net_summary_type)::tmp
    real(dp)::prec
    logical::ok
    n=size(y)
    gmax=maxval(group)
    ng=gmax
    allocate(str(n))
    str=1
    if(present(stratum)) str=stratum
    smax=maxval(str)
    ns=smax
    prec=1.0_dp; if(present(precision))prec=precision
    allocate(wr(size(times),ng,ns),we(size(times),ng,ns),dn2(size(times),ng,ns)); wr=0.0_dp; we=0.0_dp; dn2=0.0_dp
    do s=1,ns
      do g=1,ng
        k=count(group==g .and. str==s); if(k==0)cycle
        allocate(xx(k,size(x,2)),yy(k),ys(k),ev(k)); j=0
        do i=1,n
          if(group(i)==g .and. str(i)==s)then; j=j+1; xx(j,:)=x(i,:); yy(j)=y(i); ys(j)=0.0_dp; ev(j)=status(i); end if
        end do
        call netfast_summary(tab,xx,yy,ys,ev,times,prec,tmp)
        wr(:,g,s)=tmp%yisi; we(:,g,s)=tmp%dnisi-tmp%yidlisi; dn2(:,g,s)=tmp%dnisisq
        deallocate(xx,yy,ys,ev)
      end do
    end do
    allocate(totalr(size(times),ns),totale(size(times),ns),z(ng),cov(ng,ng)); z=0.0_dp; cov=0.0_dp
    do s=1,ns
      totalr(:,s)=sum(wr(:,:,s),dim=2); totale(:,s)=sum(we(:,:,s),dim=2)
      do t=1,size(times)
        if(totalr(t,s)<=0.0_dp)cycle
        do g=1,ng
          z(g)=z(g)+we(t,g,s)-wr(t,g,s)/totalr(t,s)*totale(t,s)
        end do
        do k=1,ng
          allocate(v(ng)); v=-wr(t,:,s)/totalr(t,s); v(k)=v(k)+1.0_dp
          do i=1,ng; do j=1,ng; cov(i,j)=cov(i,j)+v(i)*v(j)*dn2(t,k,s); end do; end do
          deallocate(v)
        end do
      end do
    end do
    result%df=max(0,ng-1); allocate(result%statistic_vector(result%df),result%covariance(result%df,result%df))
    if(result%df>0)then
      result%statistic_vector=z(1:result%df); result%covariance=cov(1:result%df,1:result%df)
      allocate(sol(result%df)); call solve_linear(result%covariance,result%statistic_vector,sol,ok)
      if(ok)then
        result%chisq=dot_product(result%statistic_vector,sol)
        result%p_value=chi_square_sf(result%chisq,result%df)
      end if
    end if
  end subroutine rsdiff

  subroutine nessie_expected(tab,x,group,times,expected_count,conditional_expected_time,horizon,step)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),times(:)
    integer,intent(in)::group(:)
    real(dp),allocatable,intent(out)::expected_count(:,:),conditional_expected_time(:)
    real(dp),intent(in),optional::horizon,step
    integer::ng,g,i,j,k,n
    real(dp)::hor,dt
    real(dp),allocatable::xx(:,:),tvec(:),surv(:),grid(:)
    ng=maxval(group)
    allocate(expected_count(ng,size(times)),conditional_expected_time(ng))
    expected_count=0.0_dp
    hor=100.0_dp*365.241_dp
    if(present(horizon)) hor=horizon
    dt=0.5_dp*365.241_dp
    if(present(step)) dt=step
    do g=1,ng
      k=count(group==g); if(k==0)cycle; allocate(xx(k,size(x,2))); j=0
      do i=1,size(group); if(group(i)==g)then; j=j+1; xx(j,:)=x(i,:); end if; end do
      allocate(tvec(k),surv(k))
      do j=1,size(times)
        tvec=times(j)
        surv=expected_survival_local(tab,xx,tvec)
        expected_count(g,j)=sum(surv)
      end do
      n=max(1,ceiling(hor/dt)); allocate(grid(n)); conditional_expected_time(g)=0.0_dp
      do j=1,n
        tvec=min(hor,real(j,dp)*dt)
        surv=expected_survival_local(tab,xx,tvec)
        conditional_expected_time(g)=conditional_expected_time(g) + &
          sum(surv)/real(k,dp)*dt
      end do
      deallocate(xx,tvec,surv,grid)
    end do
  contains
    function expected_survival_local(tab0,x0,t0) result(s0)
      use relsurv_ratetable, only : expected_survival
      type(ratetable_type),intent(in)::tab0
      real(dp),intent(in)::x0(:,:),t0(:)
      real(dp)::s0(size(t0))
      s0=expected_survival(tab0,x0,t0)
    end function
  end subroutine nessie_expected
end module relsurv_comparison
