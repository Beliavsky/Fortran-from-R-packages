program test_tornado_multivariate
  use mc2d, only : dp, mcnode, mc, mcdata, make_mc, tornado, tornadounc, tornado_result, &
    node_summary_each, summary_result
  implicit none

  integer, parameter :: nsv=40, nsu=30, nv=2
  real(dp) :: xa(nsv,nsu,nv), ya(nsv,nsu,nv)
  type(mcnode) :: x, y
  type(mc) :: model
  type(tornado_result) :: tv, tu
  type(summary_result),allocatable :: se(:)
  integer :: i, j, k, fails

  fails=0
  do k=1,nv
    do j=1,nsu
      do i=1,nsv
        xa(i,j,k)=(1.0_dp+0.01_dp*real(j+k,dp))*real(i,dp) &
          +real((k+1)*j,dp)
        ya(i,j,k)=2.0_dp*xa(i,j,k)+real(k,dp)
      end do
    end do
  end do

  x=mcdata(xa,type='VU',nsv=nsv,nsu=nsu,nvariates=nv,outm='each')
  y=mcdata(ya,type='VU',nsv=nsv,nsu=nsu,nvariates=nv,outm='each')
  model=make_mc([x,y],['x ','y '])

  se=node_summary_each(x)
  if(size(se)/=2)call fail('multivariate summary count',fails)

  tv=tornado(model,2)
  if(any(shape(tv%value)/=[4,2,2]))call fail('tornado multivariate shape',fails)
  if(tv%value(1,1,1)<0.999_dp)call fail('tornado multivariate correlation',fails)

  tu=tornadounc(model,2)
  if(any(shape(tu%value)/=[5,10,2]))call fail('tornadounc multivariate shape',fails)
  if(tu%value(1,1,1)<0.999_dp)call fail('tornadounc mean correlation',fails)
  if(tu%value(2,2,1)<0.999_dp)call fail('tornadounc sd correlation',fails)

  if(fails/=0)then
    print '(a,i0)','test_tornado_multivariate: FAIL ',fails
    error stop 1
  end if
  print '(a)','test_tornado_multivariate: PASS'
contains
  subroutine fail(label,nfail)
    character(len=*),intent(in)::label
    integer,intent(inout)::nfail
    print '(a)','FAIL '//trim(label)
    nfail=nfail+1
  end subroutine fail
end program test_tornado_multivariate
