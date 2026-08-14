program test_cpf
  use gpareto, only : dp, cpf_result, compute_cpf
  implicit none
  type(cpf_result)::c
  real(dp)::s1(3,4),s2(3,4),resp(3,2),ref(2)
  s1=reshape([1.0_dp,1.2_dp,0.9_dp, 2.0_dp,2.1_dp,1.8_dp, 3.0_dp,2.7_dp,3.2_dp, 4.0_dp,3.9_dp,4.1_dp],[3,4])
  s2=reshape([4.0_dp,4.1_dp,3.8_dp, 3.0_dp,2.9_dp,3.2_dp, 2.0_dp,2.2_dp,1.9_dp, 1.0_dp,1.1_dp,0.9_dp],[3,4])
  resp=reshape([1.0_dp,2.0_dp,4.0_dp, 4.0_dp,2.0_dp,1.0_dp],[3,2]);ref=[5.0_dp,5.0_dp]
  call compute_cpf(s1,s2,resp,c,ref_point=ref,n_grid=20)
  if(any(c%values<0.0_dp).or.any(c%values>1.0_dp)) error stop 'CPF probabilities invalid'
  if(c%beta_star<0.0_dp.or.c%beta_star>1.0_dp) error stop 'CPF threshold invalid'
  if(c%vd<0.0_dp) error stop 'CPF deviation invalid'
  print *, 'test_cpf PASS'
end program test_cpf
