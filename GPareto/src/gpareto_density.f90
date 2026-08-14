! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_density
  use gpareto_kinds, only : dp, i8
  use gpareto_models, only : gp_model_set, make_trend
  use gpareto_design, only : lhs_design
  use gpareto_pareto, only : nondominated_indices
  use dk_model, only : km_simulate
  implicit none
  private
  type, public :: pareto_density_result
    real(dp),allocatable :: cps(:,:),density(:)
    real(dp),allocatable :: eval_x(:,:),bandwidth(:)
  end type pareto_density_result
  public :: pareto_set_density
contains
  subroutine pareto_set_density(models,lower,upper,res,nsim,npoints,eval_points,seed)
    type(gp_model_set),intent(in)::models
    real(dp),intent(in)::lower(:),upper(:)
    type(pareto_density_result),intent(out)::res
    integer,intent(in),optional::nsim,npoints
    real(dp),intent(in),optional::eval_points(:,:)
    integer(i8),intent(in),optional::seed
    real(dp),allocatable::x(:,:),f(:,:),draw(:,:,:),tmp(:,:),cps(:,:),newcps(:,:),h(:)
    integer,allocatable::idx(:)
    integer::ns,np,d,m,i,j,k,nc
    ns=50
    if(present(nsim))ns=nsim
    np=1000
    if(present(npoints))np=npoints
    d=models%dim()
    m=models%nobj()
    call lhs_design(np,d,lower,upper,x,seed)
    allocate(draw(m,ns,np))
    do j=1,m
      call make_trend(models%model(j)%trend_kind,x,f)
      call km_simulate(models%model(j)%km,ns,x,f,tmp,conditional=.true.,nugget_sim=1.0e-8_dp)
      draw(j,:,:)=tmp
    end do
    if (allocated(tmp)) deallocate(tmp)
    allocate(cps(0,d))
    nc=0
    do i=1,ns
      allocate(tmp(np,m))
      do j=1,m
      tmp(:,j)=draw(j,i,:)
      end do
      call nondominated_indices(tmp,idx)
      allocate(newcps(nc+size(idx),d))
      if(nc>0)newcps(1:nc,:)=cps
      do k=1,size(idx)
      newcps(nc+k,:)=x(idx(k),:)
      end do
      call move_alloc(newcps,cps)
      nc=size(cps,1)
      deallocate(tmp,idx)
    end do
    res%cps=cps
    if(present(eval_points))then
    res%eval_x=eval_points
    else
    res%eval_x=x
    end if
    allocate(h(d))
    do j=1,d
    h(j)=max(stddev(cps(:,j))*real(max(nc,2),dp)**(-1.0_dp/(real(d,dp)+4.0_dp)),1.0e-8_dp)
    end do
    res%bandwidth=h
    allocate(res%density(size(res%eval_x,1)))
    call kde_diag(cps,res%eval_x,h,res%density)
  end subroutine pareto_set_density

  pure real(dp) function stddev(x) result(s)
    real(dp),intent(in)::x(:)
    real(dp)::m
    if(size(x)<2)then
    s=0.0_dp
    return
    end if
    m=sum(x)/real(size(x),dp)
    s=sqrt(sum((x-m)**2)/real(size(x)-1,dp))
  end function stddev

  subroutine kde_diag(x,q,h,dens)
    real(dp),intent(in)::x(:,:),q(:,:),h(:)
    real(dp),intent(out)::dens(:)
    real(dp)::z,norm
    integer::i,j,k,d
    d=size(x,2)
    norm=(2.0_dp*acos(-1.0_dp))**(-0.5_dp*real(d,dp))/product(h)/real(size(x,1),dp)
    do i=1,size(q,1)
    dens(i)=0.0_dp
    do j=1,size(x,1)
    z=0.0_dp
    do k=1,d
    z=z+((q(i,k)-x(j,k))/h(k))**2
    end do
    dens(i)=dens(i)+exp(-0.5_dp*z)
    end do
    dens(i)=dens(i)*norm
    end do
  end subroutine kde_diag
end module gpareto_density
