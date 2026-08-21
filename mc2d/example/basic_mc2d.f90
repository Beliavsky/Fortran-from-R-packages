program basic_mc2d
  use mc2d, only : dp, mcnode, mcdata, node_summary, summary_result, operator(+)
  implicit none
  type(mcnode) :: variability, uncertainty, total
  type(summary_result) :: s
  real(dp) :: v(100), u(50)
  integer :: i

  do i=1,size(v)
    v(i)=real(i,dp)/real(size(v),dp)
  end do
  do i=1,size(u)
    u(i)=0.2_dp*real(i,dp)/real(size(u),dp)
  end do

  variability=mcdata(v,type='V',nsv=size(v))
  uncertainty=mcdata(u,type='U',nsu=size(u))
  total=variability+uncertainty
  s=node_summary(total)

  print '(a,3(1x,i0))','total dimensions:',shape(total%value)
  print '(a,1x,a)','total type:',trim(total%type_name())
  print '(a,1x,f10.6)','median of conditional means:',s%value(1,1)
end program basic_mc2d
