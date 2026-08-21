! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_sampling
  use mc2d_kinds, only : dp, nan_dp
  use mc2d_random, only : seed_random, random_uniform_open
  use mc2d_node, only : mcnode, mc_type_0, mc_type_v, mc_type_u, mc_type_vu, ndvar, ndunc, mcdata
  implicit none
  private
  public :: mcstoc, mcprobtree_weights, mcprobtree_switch

  abstract interface
    subroutine sampler_proc(n,x)
      import dp
      integer,intent(in)::n
      real(dp),intent(out)::x(n)
    end subroutine sampler_proc
  end interface
contains
  function mcstoc(sampler,type,nsv,nsu,nvariates,outm,seed) result(node)
    procedure(sampler_proc)::sampler
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates,seed
    type(mcnode)::node
    character(len=2)::ty
    integer::nv,nu,nva,n
    real(dp),allocatable::x(:)
    ty='V';if(present(type))ty=type
    nv=ndvar();if(present(nsv))nv=nsv;nu=ndunc();if(present(nsu))nu=nsu;nva=1;if(present(nvariates))nva=nvariates
    select case(trim(ty));case('V');nu=1;case('U');nv=1;case('0');nv=1;nu=1;end select
    n=nv*nu*nva;allocate(x(n));if(present(seed))call seed_random(seed);call sampler(n,x)
    node=mcdata(x,ty,nv,nu,nva);if(present(outm))node%outm=trim(outm)
  end function mcstoc

  function mcprobtree_weights(weights,values,type,nsv,nsu,nvariates,seed) result(res)
    real(dp),intent(in)::weights(:)
    type(mcnode),intent(in)::values(:)
    character(len=*),intent(in),optional::type
    integer,intent(in),optional::nsv,nsu,nvariates,seed
    type(mcnode)::res
    integer::nv,nu,nva,ncell,i,j,k,choice
    real(dp)::u,s,c
    character(len=2)::ty
    type(mcnode),allocatable::v(:)
    if(size(weights)/=size(values).or.any(weights<0).or.sum(weights)<=0)error stop 'mcprobtree: invalid weights'
    if(present(seed))call seed_random(seed)
    ty='V'
    if(present(type))ty=type
    nv=ndvar()
    if(present(nsv))nv=nsv
    nu=ndunc()
    if(present(nsu))nu=nsu
    nva=1
    if(present(nvariates))nva=nvariates
    select case(trim(ty));case('V');nu=1;case('U');nv=1;case('0');nv=1;nu=1;end select
    allocate(v(size(values)));do i=1,size(values);v(i)=mcdata(values(i),ty,nv,nu,nva);end do
    res=mcdata(nan_dp(),ty,nv,nu,nva);s=sum(weights)
    do k=1,nva;do j=1,nu;do i=1,nv
      u=random_uniform_open()*s;c=0;choice=size(weights)
      do ncell=1,size(weights);c=c+weights(ncell);if(u<=c)then;choice=ncell;exit;end if;end do
      res%value(i,j,k)=v(choice)%value(i,j,k)
    end do;end do;end do
  end function mcprobtree_weights

  function mcprobtree_switch(switch_node,values,type,nsv,nsu,nvariates) result(res)
    type(mcnode),intent(in)::switch_node,values(:)
    character(len=*),intent(in),optional::type
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::res
    integer::nv,nu,nva,i,j,k,choice
    character(len=2)::ty
    type(mcnode),allocatable::v(:)
    ty=switch_node%type_name();if(present(type))ty=type
    nv=switch_node%nsv();if(present(nsv))nv=nsv;nu=switch_node%nsu();if(present(nsu))nu=nsu
    nva=switch_node%nvariates();if(present(nvariates))nva=nvariates
    res=mcdata(nan_dp(),ty,nv,nu,nva);allocate(v(size(values)));do i=1,size(values);v(i)=mcdata(values(i),ty,nv,nu,nva);end do
    block
      type(mcnode)::swn
      swn=mcdata(switch_node,ty,nv,nu,nva)
      do k=1,nva;do j=1,nu;do i=1,nv
        choice=nint(swn%value(i,j,k));if(choice<1.or.choice>size(values))error stop 'mcprobtree: switch out of range'
        res%value(i,j,k)=v(choice)%value(i,j,k)
      end do;end do;end do
    end block
  end function mcprobtree_switch
end module mc2d_sampling
