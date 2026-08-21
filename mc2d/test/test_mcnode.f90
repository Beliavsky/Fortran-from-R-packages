program test_mcnode
  use mc2d, only : dp, mcnode, mc, mcdata, make_mc, ndvar, ndunc, typemcnode, &
    dimmcnode, dimmc, is_mcnode, is_mc, pmin_node, pmax_node, extractvar, addvar, &
    unmc, outm_set, mcapply_reduce, mcratio, operator(+)
  implicit none
  type(mcnode) :: x0,xv,xu,xvu,r,lo,hi,e,a
  type(mc) :: model
  real(dp), allocatable :: raw(:,:,:)
  integer :: d(3), oldv, oldu, fails, i
  real(dp) :: ratio(7)

  fails=0; oldv=ndvar(); oldu=ndunc()
  d(1)=ndvar(10); d(1)=ndunc(5)
  x0=mcdata(2.0_dp,type='0')
  xv=mcdata([(real(i,dp),i=1,10)],type='V',nsv=10)
  xu=mcdata([(real(i,dp),i=1,5)],type='U',nsu=5)
  xvu=xv+xu
  if(trim(typemcnode(xvu))/='VU') call fail('V + U type',fails)
  if(any(dimmcnode(xvu)/=[10,5,1])) call fail('V + U dimensions',fails)
  call check_close(xvu%value(3,4,1),7.0_dp,'V + U value',fails)

  r=xv+x0
  if(trim(typemcnode(r))/='V') call fail('V + 0 type',fails)
  lo=pmin_node(xv,mcdata(5.0_dp,type='0'))
  hi=pmax_node(xv,mcdata(5.0_dp,type='0'))
  if(any(lo%value>5.0_dp)) call fail('pmin',fails)
  if(any(hi%value<5.0_dp)) call fail('pmax',fails)
  if(.not.is_mcnode(xv)) call fail('is_mcnode',fails)

  e=mcdata(reshape([(real(i,dp),i=1,20)],[10,1,2]),type='V',nsv=10,nvariates=2)
  r=extractvar(e,[2])
  if(r%nvariates()/=1) call fail('extractvar',fails)
  a=addvar([r,r])
  if(a%nvariates()/=2) call fail('addvar',fails)
  call outm_set(a,'mean')
  if(a%outm/='mean') call fail('outm_set',fails)
  raw=unmc(a)
  if(any(shape(raw)/=[10,1,2])) call fail('unmc',fails)

  model=make_mc([xv,xu,xvu],['xv ','xu ','sum'])
  if(.not.is_mc(model)) call fail('is_mc',fails)
  if(model%size()/=3) call fail('mc size',fails)
  d=dimmc(model)
  if(any(d/=[10,5,1])) call fail('dimmc',fails)

  r=mcapply_reduce(xvu,'unc',sum_vec)
  if(trim(typemcnode(r))/='V' .or. any(shape(r%value)/=[10,1,1])) &
    call fail('mcapply unc',fails)
  ratio=mcratio(xvu)
  if(any(ratio/=ratio)) call fail('mcratio finite',fails)

  d(1)=ndvar(oldv); d(1)=ndunc(oldu)
  if(fails/=0)then
    print '(a,i0)','test_mcnode: FAIL ',fails; error stop 1
  end if
  print '(a)','test_mcnode: PASS'
contains
  real(dp) function sum_vec(x) result(s)
    real(dp),intent(in)::x(:); s=sum(x)
  end function sum_vec
  subroutine check_close(x,y,label,fails)
    real(dp),intent(in)::x,y;character(len=*),intent(in)::label;integer,intent(inout)::fails
    if(abs(x-y)>1e-12_dp)call fail(label,fails)
  end subroutine
  subroutine fail(label,fails)
    character(len=*),intent(in)::label;integer,intent(inout)::fails
    print '(a)','FAIL '//trim(label);fails=fails+1
  end subroutine
end program test_mcnode
