program test_workflow
  use mc2d, only : dp, mcnode, mc, mcmodel, mcmodelcut, mccut_result, mcdata, make_mc, evalmcmod, &
    evalmccut_reduce, node_summary, summary_result, tornado, tornadounc, tornado_result, &
    mcprobtree_weights, mcprobtree_switch, operator(+), operator(*)
  implicit none
  type(mcmodel) :: model_def
  type(mcmodelcut) :: cut_def
  type(mccut_result) :: cut_res
  type(mc) :: model, umodel
  type(mcnode) :: x1,x2,y,u1,u2,uy,tree,switch_node
  type(summary_result) :: sm
  type(tornado_result) :: tv,tu
  real(dp) :: v1(100),v2(100),uu1(80),uu2(80)
  type(mcnode) :: choices(2)
  integer :: i,fails

  fails=0
  do i=1,100
    v1(i)=real(i,dp)
    v2(i)=sin(0.17_dp*real(i,dp))+0.01_dp*real(mod(37*i,19),dp)
  end do
  x1=mcdata(v1,type='V',nsv=100)
  x2=mcdata(v2,type='V',nsv=100)
  y=2.0_dp*x1+x2
  model=make_mc([x1,x2,y],['x1 ','x2 ','y  '])
  tv=tornado(model,3)
  if(any(tv%value < -1.0000001_dp) .or. any(tv%value > 1.0000001_dp)) &
    call fail('tornado range',fails)
  if(tv%value(1,1,1)<0.9_dp) call fail('tornado strong x1 relation',fails)

  do i=1,80
    uu1(i)=real(i,dp)
    uu2(i)=cos(0.11_dp*real(i,dp))+0.02_dp*real(mod(23*i,17),dp)
  end do
  u1=mcdata(uu1,type='U',nsu=80)
  u2=mcdata(uu2,type='U',nsu=80)
  uy=3.0_dp*u1+u2
  umodel=make_mc([u1,u2,uy],['u1 ','u2 ','uy '])
  tu=tornadounc(umodel,3)
  if(any(tu%value < -1.0000001_dp) .or. any(tu%value > 1.0000001_dp)) &
    call fail('tornadounc range',fails)

  sm=node_summary(y)
  if(size(sm%value,1)/=1 .or. size(sm%value,2)<5) &
    call fail('summary V structure',fails)

  choices(1)=mcdata(10.0_dp,type='0')
  choices(2)=mcdata(20.0_dp,type='0')
  tree=mcprobtree_weights([0.25_dp,0.75_dp],choices,type='V',nsv=200,seed=19)
  if(any(tree%value/=10.0_dp .and. tree%value/=20.0_dp)) &
    call fail('mcprobtree weighted support',fails)

  switch_node=mcdata([(real(1+mod(i,2),dp),i=1,20)],type='V',nsv=20)
  tree=mcprobtree_switch(switch_node,choices)
  if(any(tree%value/=10.0_dp .and. tree%value/=20.0_dp)) &
    call fail('mcprobtree switch support',fails)

  model_def%evaluate=>build_model
  model=evalmcmod(model_def,12,7,seed=5)
  if(model%size()/=3) call fail('evalmcmod model size',fails)
  if(any(shape(model%node(3)%value)/=[12,7,1])) &
    call fail('evalmcmod dimensions',fails)

  cut_def%evaluate_column=>build_cut_column
  cut_res=evalmccut_reduce(cut_def,10,6,reduce_cut_column,seed=7)
  if(any(shape(cut_res%value)/=[2,6]))call fail('evalmccut_reduce dimensions',fails)
  if(abs(cut_res%value(1,3)-3.55_dp)>1.0e-12_dp) &
    call fail('evalmccut_reduce mean statistic',fails)

  if(fails/=0)then
    print '(a,i0)','test_workflow: FAIL ',fails
    error stop 1
  end if
  print '(a)','test_workflow: PASS'
contains
  subroutine build_model(nsv,nsu,result)
    integer,intent(in)::nsv,nsu
    type(mc),intent(out)::result
    type(mcnode)::a,b,c
    real(dp),allocatable::av(:),bu(:)
    integer::ii
    allocate(av(nsv),bu(nsu))
    do ii=1,nsv; av(ii)=real(ii,dp)/real(nsv,dp); end do
    do ii=1,nsu; bu(ii)=real(ii,dp)/real(nsu,dp); end do
    a=mcdata(av,type='V',nsv=nsv)
    b=mcdata(bu,type='U',nsu=nsu)
    c=a+b
    result=make_mc([a,b,c],['a ','b ','c '])
  end subroutine build_model
  subroutine build_cut_column(nsv,index,result)
    integer,intent(in)::nsv,index
    type(mc),intent(out)::result
    type(mcnode)::a,b,c
    real(dp),allocatable::av(:)
    integer::ii
    allocate(av(nsv))
    do ii=1,nsv
      av(ii)=real(ii,dp)/real(nsv,dp)
    end do
    a=mcdata(av,type='V',nsv=nsv)
    b=mcdata(real(index,dp),type='0')
    c=a+b
    result=make_mc([a,b,c],['a ','b ','c '])
  end subroutine build_cut_column

  subroutine reduce_cut_column(column,statistic)
    type(mc),intent(in)::column
    real(dp),allocatable,intent(out)::statistic(:)
    allocate(statistic(2))
    statistic(1)=sum(column%node(3)%value(:,1,1))/real(column%node(3)%nsv(),dp)
    statistic(2)=maxval(column%node(3)%value(:,1,1))
  end subroutine reduce_cut_column
  subroutine fail(label,fails)
    character(len=*),intent(in)::label
    integer,intent(inout)::fails
    print '(a)','FAIL '//trim(label)
    fails=fails+1
  end subroutine fail
end program test_workflow
