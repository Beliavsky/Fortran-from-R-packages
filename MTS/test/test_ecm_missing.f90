! SPDX-License-Identifier: Artistic-2.0
program test_ecm_missing
   use mts
   use test_support
   implicit none
   integer,parameter::n=650,k=2
   real(dp)::x(n,k),u,eta,beta(k,1),sigma(k,k),piw(k,k)
   real(dp),allocatable::estimate(:),partial(:)
   logical::missing(k)
   type(vecm_model)::known,johansen
   integer::t,istat

   call set_random_seed(2222)
   x(1,:)=[0.0_dp,0.1_dp]
   do t=2,n
      u=0.20_dp*random_normal()
      eta=0.08_dp*random_normal()
      x(t,1)=x(t-1,1)+u
      x(t,2)=x(t,1)+0.55_dp*(x(t-1,2)-x(t-1,1))+eta
   end do
   beta(:,1)=[1.0_dp,-1.0_dp]
   call fit_vecm_known_beta(x,1,beta,known,include_constant=.true.)
   call assert_true(known%status==mts_success.and.known%rank==1,'known-beta VECM')
   call assert_true(all(shape(known%alpha)==[k,1]),'known-beta alpha shape')
   call fit_vecm_johansen(x,1,1,johansen,include_constant=.true.)
   call assert_true(johansen%status==mts_success.and.allocated(johansen%eigenvalues),'Johansen VECM')
   call assert_finite(johansen%eigenvalues,'finite Johansen eigenvalues')

   piw=0.0_dp;piw(1,1)=0.5_dp;piw(2,2)=0.4_dp
   sigma=reshape([0.10_dp,0.02_dp,0.02_dp,0.12_dp],[k,k])
   call estimate_missing_observation(x,piw,sigma,300,estimate,status=istat)
   call assert_true(istat==mts_success.and.size(estimate)==k,'full missing observation estimate')
   call assert_finite(estimate,'finite full missing estimate')
   missing=[.true.,.false.]
   call estimate_partial_missing(x,piw,sigma,301,missing,partial,status=istat)
   call assert_true(istat==mts_success.and.size(partial)==1,'partial missing estimate')
   call assert_finite(partial,'finite partial missing estimate')

   print '(a)','test_ecm_missing: PASS'
end program test_ecm_missing
