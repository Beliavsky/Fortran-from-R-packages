! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_cpf
  use gpareto_kinds, only : dp
  use gpareto_pareto, only : nondominated_points, hypervolume
  implicit none
  private
  type, public :: cpf_result
    real(dp), allocatable :: x(:),y(:),values(:,:),front(:,:),responses(:,:)
    real(dp), allocatable :: fun1sims(:,:),fun2sims(:,:),ve(:,:)
    real(dp) :: beta_star=-1.0_dp,vd=-1.0_dp
  end type cpf_result
  public :: compute_cpf, vorob_threshold, vorob_expectation, vorob_deviation
contains
  subroutine compute_cpf(fun1sims,fun2sims,response,res,pareto_front,f1lim,f2lim,ref_point,n_grid,compute_ve,compute_vd)
    real(dp),intent(in)::fun1sims(:,:),fun2sims(:,:),response(:,:)
    type(cpf_result),intent(out)::res
    real(dp),intent(in),optional::pareto_front(:,:),f1lim(:),f2lim(:),ref_point(:)
    integer,intent(in),optional::n_grid
    logical,intent(in),optional::compute_ve,compute_vd
    real(dp),allocatable::pf(:,:),simpts(:,:),allpts(:,:),nd(:,:)
    integer::ng,ns,np,i,j,a,b
    logical::cve,cvd,dom
    if(any(shape(fun1sims)/=shape(fun2sims))) error stop 'compute_cpf: simulation dimensions'
    if(size(response,2)/=2) error stop 'compute_cpf: two objectives required'
    ns=size(fun1sims,1)
    np=size(fun1sims,2)
    ng=100
    if(present(n_grid))ng=n_grid
    cve=.true.
    if(present(compute_ve))cve=compute_ve
    cvd=.true.
    if(present(compute_vd))cvd=compute_vd
    if(present(pareto_front))then
    pf=pareto_front
    else
    call nondominated_points(response,pf)
    end if
    if(present(f1lim))then
      res%x=f1lim
    else
      allocate(res%x(ng))
      if(present(ref_point))then
      call linspace(minval(fun1sims),ref_point(1),res%x)
      else
      call linspace(minval(fun1sims),maxval(fun1sims),res%x)
      end if
    end if
    if(present(f2lim))then
      res%y=f2lim
    else
      allocate(res%y(ng))
      if(present(ref_point))then
      call linspace(minval(fun2sims),ref_point(2),res%y)
      else
      call linspace(minval(fun2sims),maxval(fun2sims),res%y)
      end if
    end if
    allocate(res%values(size(res%x),size(res%y)))
    res%values=0.0_dp
    do a=1,size(res%x)
    do b=1,size(res%y)
      do i=1,ns
        allocate(simpts(np,2))
        simpts(:,1)=fun1sims(i,:)
        simpts(:,2)=fun2sims(i,:)
        allocate(allpts(np+size(pf,1),2))
        allpts(1:np,:)=simpts
        allpts(np+1:,:)=pf
        call nondominated_points(allpts,nd)
        dom=.false.
        do j=1,size(nd,1)
          if(all(nd(j,:)<= [res%x(a),res%y(b)]))then
          dom=.true.
          exit
          end if
        end do
        if(dom)res%values(a,b)=res%values(a,b)+1.0_dp
        deallocate(simpts,allpts,nd)
      end do
      res%values(a,b)=res%values(a,b)/real(ns,dp)
    end do
    end do
    res%front=pf
    res%responses=response
    res%fun1sims=fun1sims
    res%fun2sims=fun2sims
    if(cve)then
    res%beta_star=vorob_threshold(res)
    call vorob_expectation(res,res%ve)
    end if
    if(cvd)then
      if(.not.cve)error stop 'compute_cpf: Vorobev expectation required for deviation'
      if(present(ref_point))then
      res%vd=vorob_deviation(res,ref_point)
      else
      res%vd=vorob_deviation(res,[maxval(fun1sims),maxval(fun2sims)])
      end if
    end if
  end subroutine compute_cpf

  subroutine linspace(a,b,x)
    real(dp),intent(in)::a,b
    real(dp),intent(out)::x(:)
    integer::i,n
    n=size(x)
    if(n==1)then
    x(1)=a
    return
    end if
    do i=1,n
    x(i)=a+(b-a)*real(i-1,dp)/real(n-1,dp)
    end do
  end subroutine linspace

  real(dp) function vorob_threshold(a,precision) result(beta)
    type(cpf_result),intent(in)::a
    real(dp),intent(in),optional::precision
    real(dp)::lo,hi,pr,emu,tmp
    integer::ninside
    pr=0.001_dp
    if(present(precision))pr=precision
    lo=0.0_dp
    hi=1.0_dp
    beta=0.5_dp
    emu=sum(a%values)
    do while(hi-lo>pr)
      ninside=count(a%values>beta)
      tmp=real(ninside,dp)
      if(abs(tmp-emu)<0.5_dp)exit
      if(tmp<emu)then
      hi=beta
      else
      lo=beta
      end if
      beta=0.5_dp*(lo+hi)
    end do
  end function vorob_threshold

  subroutine vorob_expectation(a,ve)
    type(cpf_result),intent(in)::a
    real(dp),allocatable,intent(out)::ve(:,:)
    real(dp),allocatable::raw(:,:),tmp(:,:)
    integer::i,j,k,n
    allocate(raw(size(a%y),2))
    raw=huge(1.0_dp)
    k=0
    do j=size(a%y),1,-1
      n=0
      do i=1,size(a%x)
      if(a%values(i,j)<a%beta_star)n=i
      end do
      if(n>0.and.n<size(a%x))then
      k=k+1
      raw(k,:)=[a%x(n),a%y(j)]
      end if
    end do
    if(k==0)then
    allocate(ve(0,2))
    return
    end if
    call nondominated_points(raw(1:k,:),tmp)
    ve=tmp
  end subroutine vorob_expectation

  real(dp) function vorob_deviation(a,ref) result(vd)
    type(cpf_result),intent(in)::a
    real(dp),intent(in)::ref(2)
    real(dp),allocatable::s(:,:),u(:,:)
    real(dp)::h1,h2,h12
    integer::i,ns,np
    ns=size(a%fun1sims,1)
    np=size(a%fun1sims,2)
    h2=hypervolume(a%ve,ref)
    vd=0.0_dp
    do i=1,ns
      allocate(s(np,2))
      s(:,1)=a%fun1sims(i,:)
      s(:,2)=a%fun2sims(i,:)
      h1=hypervolume(s,ref)
      allocate(u(np+size(a%ve,1),2))
      u(1:np,:)=s
      u(np+1:,:)=a%ve
      h12=hypervolume(u,ref)
      vd=vd+2*h12-h1-h2
      deallocate(s,u)
    end do
    vd=vd/real(ns,dp)
  end function vorob_deviation
end module gpareto_cpf
