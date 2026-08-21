module relsurv_nonparametric
  use relsurv_kinds, only : dp
  use relsurv_ratetable, only : ratetable_type, net_summary_type, netfast_summary, netwei_summary
  implicit none
  private
  integer, parameter, public :: method_pohar_perme=1, method_ederer2=2, method_hakulinen=3, method_ederer1=4
  integer, parameter, public :: transform_km=1, transform_fh=2
  integer, parameter, public :: conf_plain=1, conf_log=2, conf_loglog=3
  type, public :: rs_surv_result
    real(dp), allocatable :: time(:), n_risk(:), n_event(:), n_censor(:)
    real(dp), allocatable :: surv(:), std_err(:), lower(:), upper(:)
  end type rs_surv_result
  public :: rs_surv, cumulative_incidence_competing
contains
  subroutine rs_surv(tab,x,y,status,times,result,method,transform,ystart,precision,potential_time,conf_type,zcrit)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::x(:,:),y(:),times(:)
    integer,intent(in)::status(:)
    type(rs_surv_result),intent(out)::result
    integer,intent(in),optional::method,transform,conf_type
    real(dp),intent(in),optional::ystart(:),precision,potential_time(:),zcrit
    type(net_summary_type)::s,s2
    real(dp),allocatable::start(:),pot(:),haz(:),popsur(:)
    integer,allocatable::zero_status(:)
    real(dp)::prec,z
    integer::m,tr,ct,j,n
    n=size(y); m=method_pohar_perme; if(present(method))m=method
    tr=transform_km; if(present(transform))tr=transform
    ct=conf_log; if(present(conf_type))ct=conf_type
    prec=1.0_dp; if(present(precision))prec=precision
    z=1.959963984540054_dp; if(present(zcrit))z=zcrit
    allocate(start(n)); start=0.0_dp; if(present(ystart))start=ystart
    allocate(result%time(size(times)),result%n_risk(size(times)),result%n_event(size(times)), &
      result%n_censor(size(times)),result%surv(size(times)),result%std_err(size(times)), &
      result%lower(size(times)),result%upper(size(times)),haz(size(times)),popsur(size(times)))
    result%time=times; haz=0.0_dp; popsur=1.0_dp
    if(m<=method_ederer2) then
      call netfast_summary(tab,x,y,start,status,times,prec,s)
    else
      call netwei_summary(tab,x,y,status,times,s)
    end if
    result%n_risk=s%yi; result%n_event=s%dni; result%n_censor=s%dci
    select case(m)
    case(method_pohar_perme)
      do j=1,size(times)
        if(s%yisi(j)>0.0_dp) then
          haz(j)=s%dnisi(j)/s%yisi(j)-s%yidlisiw(j)
          result%std_err(j)=sqrt(sum(s%dnisisq(1:j)/max(s%yisi(1:j)**2,tiny(1.0_dp))))
        else
          haz(j)=0.0_dp; result%std_err(j)=0.0_dp
        end if
      end do
    case(method_ederer2)
      do j=1,size(times)
        if(s%yi(j)>0.0_dp) haz(j)=(s%dni(j)-s%yidli(j))/s%yi(j)
        result%std_err(j)=sqrt(sum(s%dni(1:j)/max(s%yi(1:j)**2,tiny(1.0_dp))))
      end do
    case(method_hakulinen)
      if(.not.present(potential_time)) error stop 'rs_surv: potential_time required for Hakulinen method'
      allocate(pot(n),zero_status(n)); pot=potential_time; zero_status=0
      call netwei_summary(tab,x,pot,zero_status,times,s2)
      popsur=1.0_dp
      do j=1,size(times)
        if(s2%yisis(j)>0.0_dp) then
          if(j==1) then
            popsur(j)=exp(-s2%yisidli(j)/s2%yisis(j))
          else
            popsur(j)=popsur(j-1)*exp(-s2%yisidli(j)/s2%yisis(j))
          end if
        else if(j>1) then
          popsur(j)=popsur(j-1)
        end if
        if(s%yi(j)>0.0_dp) haz(j)=s%dni(j)/s%yi(j)
        result%std_err(j)=sqrt(sum(s%dni(1:j)/max(s%yi(1:j)**2,tiny(1.0_dp))))
      end do
    case(method_ederer1)
      popsur=s%sis/real(n,dp)
      do j=1,size(times)
        if(s%yi(j)>0.0_dp) haz(j)=s%dni(j)/s%yi(j)
        result%std_err(j)=sqrt(sum(s%dni(1:j)/max(s%yi(1:j)**2,tiny(1.0_dp))))
      end do
    case default
      error stop 'rs_surv: unknown method'
    end select
    if(tr==transform_fh) then
      do j=1,size(times)
        result%surv(j)=exp(-sum(haz(1:j)))
      end do
    else
      result%surv=1.0_dp
      do j=1,size(times)
        if(j==1) then
          result%surv(j)=1.0_dp-haz(j)
        else
          result%surv(j)=result%surv(j-1)*(1.0_dp-haz(j))
        end if
      end do
    end if
    if(m>method_ederer2) result%surv=result%surv/max(popsur,tiny(1.0_dp))
    call make_ci(result%surv,result%std_err,ct,z,result%lower,result%upper)
  end subroutine rs_surv

  subroutine make_ci(s,se,ct,z,lo,hi)
    real(dp),intent(in)::s(:),se(:),z
    integer,intent(in)::ct
    real(dp),intent(out)::lo(:),hi(:)
    integer::i
    do i=1,size(s)
      select case(ct)
      case(conf_plain)
        lo(i)=s(i)*(1.0_dp-z*se(i)); hi(i)=s(i)*(1.0_dp+z*se(i))
      case(conf_loglog)
        if(s(i)>0.0_dp .and. s(i)<1.0_dp) then
          lo(i)=exp(-exp(log(-log(s(i)))-z*se(i)/log(s(i))))
          hi(i)=exp(-exp(log(-log(s(i)))+z*se(i)/log(s(i))))
        else
          lo(i)=s(i); hi(i)=s(i)
        end if
      case default
        if(s(i)>0.0_dp) then
          lo(i)=exp(log(s(i))-z*se(i)); hi(i)=exp(log(s(i))+z*se(i))
        else
          lo(i)=0.0_dp; hi(i)=0.0_dp
        end if
      end select
    end do
  end subroutine make_ci

  subroutine cumulative_incidence_competing(time,nrisk,event_cause,cause,cif,surv)
    real(dp),intent(in)::time(:),nrisk(:)
    integer,intent(in)::event_cause(:,:),cause
    real(dp),intent(out)::cif(:),surv(:)
    real(dp)::s0,dh_all,dh_c,cif0
    integer::j
    s0=1.0_dp; cif0=0.0_dp; cif=0.0_dp; surv=1.0_dp
    do j=1,size(time)
      if(nrisk(j)>0.0_dp) then
        dh_all=real(count(event_cause(:,j)>0),dp)/nrisk(j)
        dh_c=real(count(event_cause(:,j)==cause),dp)/nrisk(j)
      else
        dh_all=0.0_dp; dh_c=0.0_dp
      end if
      cif0=cif0+s0*dh_c
      cif(j)=cif0
      s0=s0*(1.0_dp-dh_all)
      surv(j)=s0
    end do
  end subroutine cumulative_incidence_competing
end module relsurv_nonparametric
