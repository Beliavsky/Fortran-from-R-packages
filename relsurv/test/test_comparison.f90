program test_comparison
  use relsurv
  implicit none
  type(ratetable_type) :: tab
  type(cmp_rel_result) :: cr
  type(rsdiff_result) :: rd
  integer :: dims(1),factor(1),ncuts(1)
  integer :: status2(2),status4(4),group4(4),groupn(4)
  real(dp) :: cuts(1,1),rates(1),x2(2,1),y2(2),times2(2)
  real(dp) :: x4(4,1),y4(4),times4(2),xn(4,1),tn(2),lambda
  real(dp), allocatable :: ec(:,:),cet(:)

  dims=1; factor=1; ncuts=0; cuts=0.0_dp; rates=0.0_dp
  tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  x2=1.0_dp; y2=[1.0_dp,2.0_dp]; status2=[1,0]; times2=y2
  call cmp_rel(tab,x2,y2,status2,times2,cr,scale=1.0_dp)
  call assert_close(cr%disease(1),0.5_dp,1.0e-12_dp,'cmp disease 1')
  call assert_close(cr%disease(2),0.5_dp,1.0e-12_dp,'cmp disease 2')
  call assert_close(maxval(abs(cr%population)),0.0_dp,1.0e-12_dp,'cmp population')

  x4=1.0_dp; y4=[1.0_dp,2.0_dp,1.0_dp,2.0_dp]
  status4=[1,0,1,0]; group4=[1,1,2,2]; times4=[1.0_dp,2.0_dp]
  call rsdiff(tab,x4,y4,status4,group4,times4,rd)
  call assert_close(rd%chisq,0.0_dp,1.0e-12_dp,'rsdiff balanced')

  lambda=0.01_dp; rates=lambda
  tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  xn=1.0_dp; groupn=[1,1,2,2]; tn=[1.0_dp,10.0_dp]
  call nessie_expected(tab,xn,groupn,tn,ec,cet,horizon=10.0_dp,step=1.0_dp)
  call assert_close(ec(1,1),2.0_dp*exp(-lambda),1.0e-12_dp,'nessie g1 t1')
  call assert_close(ec(2,2),2.0_dp*exp(-10.0_dp*lambda),1.0e-12_dp,'nessie g2 t2')

  print *, 'test_comparison: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol) then
      print *, 'FAIL ',trim(msg),a,b
      error stop 1
    end if
  end subroutine assert_close
end program test_comparison
