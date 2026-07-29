! SPDX-License-Identifier: GPL-3.0-or-later
program test_copula_risk
  use nvmix
  implicit none
  real(dp) :: scale(2,2),v,loc1(1),scale1(1,1),dfs(1)
  type(nvmix_model) :: model
  scale=reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
  v=dstudent_copula([0.2_dp,0.7_dp],7.0_dp,scale)
  call assert_close(v,0.7711473286664019_dp,2.0e-11_dp,'t copula density')
  loc1=0.0_dp; scale1(1,1)=1.0_dp
  model=make_nvmix_model(loc1,scale1,mix_inverse_gamma,7.0_dp)
  call assert_close(var_nvmix(0.99_dp,model),2.9979515668685277_dp,2.0e-10_dp,'Student VaR')
  call assert_close(es_nvmix(0.99_dp,model),3.7699267861721952_dp,2.0e-10_dp,'Student ES')
  dfs=[7.0_dp]
  call assert_close(lambda_gstudent(dfs,0.4_dp),0.10122151579386639_dp,2.0e-12_dp,'tail dependence')
  call assert_close(kendall_nvmix(0.5_dp),1.0_dp/3.0_dp,2.0e-15_dp,'Kendall tau')
  print '(a)','test_copula_risk: PASS'
contains
  subroutine assert_close(a,b,tol,label)
    real(dp), intent(in) :: a,b,tol
    character(*), intent(in) :: label
    if(abs(a-b)>tol*max(1.0_dp,abs(b)))then
      write(*,'(a,3es24.15)')trim(label)//' mismatch: ',a,b,abs(a-b); error stop 1
    end if
  end subroutine
end program
