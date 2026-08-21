program test_derivatives
  use gb2, only : dp,logf_gb2,dlogf_gb2,d2logf_gb2,info_gb2
  implicit none
  real(dp)::par(4),g(4),h(4,4),gn(4),hn(4,4),pp(4),pm(4),gp(4),gm(4),step,info(4,4)
  integer::j,fails
  par=[2.3_dp,4.2_dp,1.7_dp,3.4_dp]
  fails=0
  call dlogf_gb2(3.1_dp,par(1),par(2),par(3),par(4),g)
  call d2logf_gb2(3.1_dp,par(1),par(2),par(3),par(4),h)
  do j=1,4
    step=1e-6_dp*max(1._dp,abs(par(j)))
    pp=par
    pm=par
    pp(j)=pp(j)+step
    pm(j)=pm(j)-step
    gn(j)=(logf_gb2(3.1_dp,pp(1),pp(2),pp(3),pp(4))-logf_gb2(3.1_dp,pm(1),pm(2),pm(3),pm(4)))/(2*step)
    call dlogf_gb2(3.1_dp,pp(1),pp(2),pp(3),pp(4),gp)
    call dlogf_gb2(3.1_dp,pm(1),pm(2),pm(3),pm(4),gm)
    hn(:,j)=(gp-gm)/(2*step)
  end do
  if(maxval(abs(g-gn))>2e-7_dp) then
  print *,'score error',maxval(abs(g-gn))
  fails=fails+1
  end if
  if(maxval(abs(h-hn))>2e-6_dp) then
  print *,'hessian error',maxval(abs(h-hn))
  fails=fails+1
  end if
  if(maxval(abs(h-transpose(h)))>1e-12_dp) fails=fails+1
  call info_gb2(par(1),par(2),par(3),par(4),info)
  if(maxval(abs(info-transpose(info)))>1e-12_dp) fails=fails+1
  if(fails>0) error stop 1
  print '(a)','test_derivatives: PASS'
end program
