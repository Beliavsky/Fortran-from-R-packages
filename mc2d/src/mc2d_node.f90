! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_node
  use mc2d_kinds, only : dp, nan_dp
  implicit none
  private

  integer, parameter, public :: mc_type_0=0, mc_type_v=1, mc_type_u=2, mc_type_vu=3
  integer, save :: default_nsv=1001, default_nsu=101

  type, public :: mcnode
    real(dp), allocatable :: value(:,:,:)
    integer :: node_type=mc_type_0
    character(len=:), allocatable :: outm
  contains
    procedure :: nsv => node_nsv
    procedure :: nsu => node_nsu
    procedure :: nvariates => node_nvariates
    procedure :: type_name => node_type_name
  end type mcnode

  type, public :: mc
    type(mcnode), allocatable :: node(:)
    character(len=64), allocatable :: name(:)
  contains
    procedure :: size => mc_size
  end type mc

  public :: ndvar, ndunc, mcdata, mcdatanocontrol, make_mc
  public :: dimmcnode, dimmc, typemcnode, is_mcnode, is_mc
  public :: pmin_node, pmax_node, extractvar, addvar, unmc, outm_set
  public :: operator(+),operator(-),operator(*),operator(/),operator(**)

  interface mcdata
    module procedure mcdata_scalar,mcdata_vector,mcdata_matrix,mcdata_array,mcdata_from_node
  end interface
  interface operator(+)
    module procedure node_add_node,node_add_scalar,scalar_add_node
  end interface
  interface operator(-)
    module procedure node_sub_node,node_sub_scalar,scalar_sub_node,node_neg
  end interface
  interface operator(*)
    module procedure node_mul_node,node_mul_scalar,scalar_mul_node
  end interface
  interface operator(/)
    module procedure node_div_node,node_div_scalar,scalar_div_node
  end interface
  interface operator(**)
    module procedure node_pow_scalar
  end interface
contains

  integer function ndvar(n) result(v)
    integer,intent(in),optional::n
    if(present(n))then
      if(n<=0)error stop 'ndvar: n must be positive'
      default_nsv=n
    end if
    v=default_nsv
  end function ndvar
  integer function ndunc(n) result(v)
    integer,intent(in),optional::n
    if(present(n))then
      if(n<=0)error stop 'ndunc: n must be positive'
      default_nsu=n
    end if
    v=default_nsu
  end function ndunc

  integer function node_nsv(self) result(n)
    class(mcnode),intent(in)::self
    if(allocated(self%value))then;n=size(self%value,1);else;n=0;end if
  end function
  integer function node_nsu(self) result(n)
    class(mcnode),intent(in)::self
    if(allocated(self%value))then;n=size(self%value,2);else;n=0;end if
  end function
  integer function node_nvariates(self) result(n)
    class(mcnode),intent(in)::self
    if(allocated(self%value))then;n=size(self%value,3);else;n=0;end if
  end function
  character(len=2) function node_type_name(self) result(s)
    class(mcnode),intent(in)::self
    select case(self%node_type)
    case(mc_type_0);s='0 '
    case(mc_type_v);s='V '
    case(mc_type_u);s='U '
    case(mc_type_vu);s='VU'
    case default;s='? '
    end select
  end function
  integer function mc_size(self) result(n)
    class(mc),intent(in)::self
    if(allocated(self%node))then;n=size(self%node);else;n=0;end if
  end function

  subroutine dims_for_type(type_name,nsv,nsu,nvariates,d)
    character(len=*),intent(in)::type_name
    integer,intent(in)::nsv,nsu,nvariates
    integer,intent(out)::d(3)
    select case(trim(type_name))
    case('V');d=[nsv,1,nvariates]
    case('U');d=[1,nsu,nvariates]
    case('VU');d=[nsv,nsu,nvariates]
    case('0');d=[1,1,nvariates]
    case default;error stop 'mcdata: invalid node type'
    end select
  end subroutine
  integer function type_code(type_name) result(t)
    character(len=*),intent(in)::type_name
    select case(trim(type_name));case('0');t=mc_type_0;case('V');t=mc_type_v
    case('U');t=mc_type_u;case('VU');t=mc_type_vu;case default;t=-1;end select
  end function

  function mcdata_scalar(data,type,nsv,nsu,nvariates,outm) result(x)
    real(dp),intent(in)::data
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::x
    character(len=2)::ty;integer::nv,nu,nva,d(3)
    ty='V';if(present(type))ty=type;nv=ndvar();if(present(nsv))nv=nsv
    nu=ndunc();if(present(nsu))nu=nsu;nva=1;if(present(nvariates))nva=nvariates
    call dims_for_type(ty,nv,nu,nva,d);allocate(x%value(d(1),d(2),d(3)));x%value=data
    x%node_type=type_code(ty);x%outm='each';if(present(outm))x%outm=trim(outm)
  end function

  function mcdata_vector(data,type,nsv,nsu,nvariates,outm) result(x)
    real(dp),intent(in)::data(:)
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::x
    real(dp),allocatable::flat(:)
    character(len=2)::ty;integer::nv,nu,nva,d(3),ncell,i
    ty='V';if(present(type))ty=type;nv=ndvar();if(present(nsv))nv=nsv
    nu=ndunc();if(present(nsu))nu=nsu;nva=1;if(present(nvariates))nva=nvariates
    call dims_for_type(ty,nv,nu,nva,d);ncell=product(d)
    if(size(data)/=1 .and. size(data)/=d(1)*d(2) .and. size(data)/=ncell) &
      error stop 'mcdata: vector size incompatible with node dimensions'
    allocate(x%value(d(1),d(2),d(3)),flat(ncell))
    do i=1,ncell
      flat(i)=data(1+mod(i-1,size(data)))
    end do
    x%value=reshape(flat,d)
    x%node_type=type_code(ty);x%outm='each';if(present(outm))x%outm=trim(outm)
  end function

  function mcdata_matrix(data,type,nsv,nsu,nvariates,outm) result(x)
    real(dp),intent(in)::data(:,:)
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::x
    real(dp),allocatable::flat(:)
    allocate(flat(size(data)));flat=reshape(data,[size(data)])
    x=mcdata_vector(flat,type,nsv,nsu,nvariates,outm)
  end function

  function mcdata_array(data,type,nsv,nsu,nvariates,outm) result(x)
    real(dp),intent(in)::data(:,:,:)
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::x
    character(len=2)::ty;integer::nv,nu,nva,d(3),i,j,k
    ty='V';if(present(type))ty=type;nv=ndvar();if(present(nsv))nv=nsv
    nu=ndunc();if(present(nsu))nu=nsu;nva=size(data,3);if(present(nvariates))nva=nvariates
    call dims_for_type(ty,nv,nu,nva,d)
    if(.not.(size(data,1)==d(1) .or. size(data,1)==1))error stop 'mcdata: incompatible first dimension'
    if(.not.(size(data,2)==d(2) .or. size(data,2)==1))error stop 'mcdata: incompatible second dimension'
    if(.not.(size(data,3)==d(3) .or. size(data,3)==1))error stop 'mcdata: incompatible variates dimension'
    allocate(x%value(d(1),d(2),d(3)))
    do k=1,d(3);do j=1,d(2);do i=1,d(1)
      x%value(i,j,k)=data(1+mod(i-1,size(data,1)),1+mod(j-1,size(data,2)),1+mod(k-1,size(data,3)))
    end do;end do;end do
    x%node_type=type_code(ty);x%outm='each';if(present(outm))x%outm=trim(outm)
  end function

  function mcdata_from_node(data,type,nsv,nsu,nvariates,outm) result(x)
    type(mcnode),intent(in)::data
    character(len=*),intent(in),optional::type,outm
    integer,intent(in),optional::nsv,nsu,nvariates
    type(mcnode)::x
    character(len=2)::ty;integer::nv,nu,nva
    ty=data%type_name();if(present(type))ty=type
    nv=max(1,data%nsv());if(present(nsv))nv=nsv;nu=max(1,data%nsu());if(present(nsu))nu=nsu
    nva=max(1,data%nvariates());if(present(nvariates))nva=nvariates
    x=mcdata_array(data%value,ty,nv,nu,nva)
    if(allocated(data%outm))x%outm=data%outm;if(present(outm))x%outm=trim(outm)
  end function

  function mcdatanocontrol(data,type,nsv,nsu,nvariates,outm) result(x)
    real(dp),intent(in)::data(:);character(len=*),intent(in)::type
    integer,intent(in)::nsv,nsu,nvariates;character(len=*),intent(in),optional::outm
    type(mcnode)::x;integer::i,n
    real(dp),allocatable::flat(:)
    allocate(x%value(nsv,nsu,nvariates));n=product(shape(x%value));allocate(flat(n))
    do i=1,n
      flat(i)=data(1+mod(i-1,size(data)))
    end do
    x%value=reshape(flat,[nsv,nsu,nvariates])
    x%node_type=type_code(type);x%outm='each';if(present(outm))x%outm=trim(outm)
  end function

  function make_mc(nodes,names) result(m)
    type(mcnode),intent(in)::nodes(:);character(len=*),intent(in),optional::names(:)
    type(mc)::m;integer::i,nsvm,nsum
    nsvm=maxval([(nodes(i)%nsv(),i=1,size(nodes))]);nsum=maxval([(nodes(i)%nsu(),i=1,size(nodes))])
    do i=1,size(nodes)
      if(nodes(i)%nsv()/=1 .and. nodes(i)%nsv()/=nsvm)error stop 'make_mc: inconsistent variability dimensions'
      if(nodes(i)%nsu()/=1 .and. nodes(i)%nsu()/=nsum)error stop 'make_mc: inconsistent uncertainty dimensions'
    end do
    allocate(m%node(size(nodes)),m%name(size(nodes)));m%node=nodes
    do i=1,size(nodes);write(m%name(i),'("node",i0)')i;end do
    if(present(names))then
      if(size(names)/=size(nodes))error stop 'make_mc: names size mismatch'
      do i=1,size(nodes);m%name(i)=names(i);end do
    end if
  end function

  function dimmcnode(x) result(d)
    type(mcnode),intent(in)::x;integer::d(3);d=shape(x%value)
  end function
  function dimmc(x) result(d)
    type(mc),intent(in)::x;integer::d(3),i;d=0
    do i=1,x%size();d=max(d,shape(x%node(i)%value));end do
  end function
  character(len=2) function typemcnode(x) result(s)
    type(mcnode),intent(in)::x;s=x%type_name()
  end function
  logical function is_mcnode(x) result(ok)
    type(mcnode),intent(in)::x
    ok=allocated(x%value) .and. x%node_type>=mc_type_0 .and. x%node_type<=mc_type_vu
    if(.not.ok)return
    select case(x%node_type);case(mc_type_0);ok=x%nsv()==1.and.x%nsu()==1
    case(mc_type_v);ok=x%nsu()==1;case(mc_type_u);ok=x%nsv()==1;case(mc_type_vu);ok=.true.;end select
  end function
  logical function is_mc(x) result(ok)
    type(mc),intent(in)::x;integer::i,d(3)
    ok=allocated(x%node);if(.not.ok)return;d=dimmc(x)
    do i=1,size(x%node)
      if(.not.is_mcnode(x%node(i)))then;ok=.false.;return;end if
      if(x%node(i)%nsv()/=1 .and. x%node(i)%nsv()/=d(1))ok=.false.
      if(x%node(i)%nsu()/=1 .and. x%node(i)%nsu()/=d(2))ok=.false.
    end do
  end function

  subroutine broadcast_pair(a,b,av,bv,d)
    type(mcnode),intent(in)::a,b;real(dp),allocatable,intent(out)::av(:,:,:),bv(:,:,:);integer,intent(out)::d(3)
    integer::da(3),db(3),i,j,k
    da=shape(a%value);db=shape(b%value)
    if(any((da/=db).and.(da/=1).and.(db/=1)))error stop 'mcnode operation: incompatible dimensions'
    d=max(da,db);allocate(av(d(1),d(2),d(3)),bv(d(1),d(2),d(3)))
    do k=1,d(3);do j=1,d(2);do i=1,d(1)
      av(i,j,k)=a%value(1+mod(i-1,da(1)),1+mod(j-1,da(2)),1+mod(k-1,da(3)))
      bv(i,j,k)=b%value(1+mod(i-1,db(1)),1+mod(j-1,db(2)),1+mod(k-1,db(3)))
    end do;end do;end do
  end subroutine
  subroutine finish_binary(a,b,d,r)
    type(mcnode),intent(in)::a,b;integer,intent(in)::d(3);type(mcnode),intent(inout)::r
    if(d(1)==1.and.d(2)==1)then;r%node_type=mc_type_0
    else if(d(2)==1)then;r%node_type=mc_type_v
    else if(d(1)==1)then;r%node_type=mc_type_u
    else;r%node_type=mc_type_vu;end if
    if(allocated(a%outm))then;r%outm=a%outm;else if(allocated(b%outm))then;r%outm=b%outm;else;r%outm='each';end if
  end subroutine

  function node_add_node(a,b) result(r)
  type(mcnode),intent(in)::a,b
  type(mcnode)::r
  real(dp),allocatable::av(:,:,:),bv(:,:,:)
  integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=av+bv;call finish_binary(a,b,d,r);end function
  function node_sub_node(a,b) result(r)
  type(mcnode),intent(in)::a,b
  type(mcnode)::r
  real(dp),allocatable::av(:,:,:),bv(:,:,:)
  integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=av-bv;call finish_binary(a,b,d,r);end function
  function node_mul_node(a,b) result(r)
  type(mcnode),intent(in)::a,b
  type(mcnode)::r
  real(dp),allocatable::av(:,:,:),bv(:,:,:)
  integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=av*bv;call finish_binary(a,b,d,r);end function
  function node_div_node(a,b) result(r)
  type(mcnode),intent(in)::a,b
  type(mcnode)::r
  real(dp),allocatable::av(:,:,:),bv(:,:,:)
  integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=av/bv;call finish_binary(a,b,d,r);end function
  function node_add_scalar(a,b) result(r)
  type(mcnode),intent(in)::a
  real(dp),intent(in)::b
  type(mcnode)::r
  r=a
  r%value=a%value+b
  end function
  function scalar_add_node(a,b) result(r)
  real(dp),intent(in)::a
  type(mcnode),intent(in)::b
  type(mcnode)::r
  r=b
  r%value=a+b%value
  end function
  function node_sub_scalar(a,b) result(r)
  type(mcnode),intent(in)::a
  real(dp),intent(in)::b
  type(mcnode)::r
  r=a
  r%value=a%value-b
  end function
  function scalar_sub_node(a,b) result(r)
  real(dp),intent(in)::a
  type(mcnode),intent(in)::b
  type(mcnode)::r
  r=b
  r%value=a-b%value
  end function
  function node_mul_scalar(a,b) result(r)
  type(mcnode),intent(in)::a
  real(dp),intent(in)::b
  type(mcnode)::r
  r=a
  r%value=a%value*b
  end function
  function scalar_mul_node(a,b) result(r)
  real(dp),intent(in)::a
  type(mcnode),intent(in)::b
  type(mcnode)::r
  r=b
  r%value=a*b%value
  end function
  function node_div_scalar(a,b) result(r)
  type(mcnode),intent(in)::a
  real(dp),intent(in)::b
  type(mcnode)::r
  r=a
  r%value=a%value/b
  end function
  function scalar_div_node(a,b) result(r)
  real(dp),intent(in)::a
  type(mcnode),intent(in)::b
  type(mcnode)::r
  r=b
  r%value=a/b%value
  end function
  function node_pow_scalar(a,b) result(r)
  type(mcnode),intent(in)::a
  real(dp),intent(in)::b
  type(mcnode)::r
  r=a
  r%value=a%value**b
  end function
  function node_neg(a) result(r);type(mcnode),intent(in)::a;type(mcnode)::r;r=a;r%value=-a%value;end function

  function pmin_node(a,b) result(r)
    type(mcnode),intent(in)::a,b;type(mcnode)::r;real(dp),allocatable::av(:,:,:),bv(:,:,:);integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=min(av,bv);call finish_binary(a,b,d,r)
  end function
  function pmax_node(a,b) result(r)
    type(mcnode),intent(in)::a,b;type(mcnode)::r;real(dp),allocatable::av(:,:,:),bv(:,:,:);integer::d(3)
    call broadcast_pair(a,b,av,bv,d);allocate(r%value(d(1),d(2),d(3)));r%value=max(av,bv);call finish_binary(a,b,d,r)
  end function

  function extractvar(x,which) result(r)
    type(mcnode),intent(in)::x;integer,intent(in)::which(:);type(mcnode)::r;integer::k
    if(any(which<1).or.any(which>x%nvariates()))error stop 'extractvar: invalid variate index'
    allocate(r%value(x%nsv(),x%nsu(),size(which)))
    do k=1,size(which);r%value(:,:,k)=x%value(:,:,which(k));end do
    r%node_type=x%node_type;if(allocated(x%outm))r%outm=x%outm
  end function

  function addvar(nodes) result(r)
    type(mcnode),intent(in)::nodes(:);type(mcnode)::r;integer::i,k0,k1,nva
    if(size(nodes)==0)error stop 'addvar: empty input';nva=sum([(nodes(i)%nvariates(),i=1,size(nodes))])
    do i=2,size(nodes)
      if(nodes(i)%nsv()/=nodes(1)%nsv().or.nodes(i)%nsu()/=nodes(1)%nsu().or.nodes(i)%node_type/=nodes(1)%node_type) &
        error stop 'addvar: nodes must have same first two dimensions and type'
    end do
    allocate(r%value(nodes(1)%nsv(),nodes(1)%nsu(),nva));k0=1
    do i=1,size(nodes);k1=k0+nodes(i)%nvariates()-1;r%value(:,:,k0:k1)=nodes(i)%value;k0=k1+1;end do
    r%node_type=nodes(1)%node_type;if(allocated(nodes(1)%outm))r%outm=nodes(1)%outm
  end function

  function unmc(x) result(a)
    type(mcnode),intent(in)::x;real(dp),allocatable::a(:,:,:);a=x%value
  end function
  subroutine outm_set(x,value)
    type(mcnode),intent(inout)::x;character(len=*),intent(in)::value;x%outm=trim(value)
  end subroutine
end module mc2d_node
