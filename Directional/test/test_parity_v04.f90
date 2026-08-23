program test_parity_v04
   use directional
   implicit none
   real(dp)::x(8,3),u(8),nrm
   integer::i,fail
   type(vector_mle_result)::a,b,c,d
   type(test_result)::t
   fail=0
   do i=1,8
      x(i,:)=[1.0_dp,0.12_dp*cos(real(i,dp)),0.12_dp*sin(real(i,dp))]
      nrm=sqrt(sum(x(i,:)**2));x(i,:)=x(i,:)/nrm
      u(i)=10.0_dp+2.0_dp*real(i-1,dp)
   end do
   a=iag_mle(x,maxit=200);if(.not.(a%loglik>-huge(1.0_dp).and.size(a%mu)==3))fail=fail+1
   b=sipc_mle(x,maxit=200);if(.not.(b%loglik>-huge(1.0_dp).and.b%gamma>=0.0_dp))fail=fail+1
   c=cipc_mle(u);if(.not.(c%gamma>=0.0_dp.and.c%angle>=0.0_dp.and.c%angle<360.0_dp))fail=fail+1
   d=gcpc_mle(u);if(.not.(d%mu(3)>0.0_dp.and.d%gamma>=0.0_dp))fail=fail+1
   t=embed_perm(x(1:4,:),x(5:8,:),19);if(t%p_value<=0.0_dp.or.t%p_value>1.0_dp)fail=fail+1
   t=hcf_perm(x(1:4,:),x(5:8,:),19);if(t%p_value<=0.0_dp.or.t%p_value>1.0_dp)fail=fail+1
   t=het_perm(x(1:4,:),x(5:8,:),9);if(t%p_value<=0.0_dp.or.t%p_value>1.0_dp)fail=fail+1
   if(fail==0)then
      print '(a)','test_parity_v04: PASS'
   else
      print '(a,i0)','test_parity_v04: FAIL ',fail;error stop 1
   end if
end program
