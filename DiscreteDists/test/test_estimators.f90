program test_estimators
   use discretedists
   implicit none
   real(dp)::y(12),par(2),mu
   integer::st,fails
   type(discrete_family_t)::f
   y=[0._dp,0._dp,1._dp,1._dp,1._dp,2._dp,2._dp,3._dp,3._dp,4._dp,5._dp,6._dp]
   fails=0
   par=estim_mu_sigma_compo(y)
   call chk(all(par>0.0_dp),'COMPO initial')
   mu=estim_mu_dbh(y,st);call chk(mu>0.0_dp.and.mu<1.0_dp,'DBH estimate')
   mu=estim_mu_dld(y,st);call chk(mu>0.0_dp,'DLD estimate')
   par=estim_mu_sigma_ggeo(y,st)
   call chk(par(1)>0.0_dp.and.par(1)<1.0_dp.and.par(2)>0.0_dp,'GGEO estimate')
   par=estim_mu_sigma_dgeii(y,st)
   call chk(par(1)>0.0_dp.and.par(1)<1.0_dp.and.par(2)>0.0_dp,'DGEII estimate')
   f=dperks();call f%initial(y,mu,par(2),st)
   call chk(abs(mu-1.0_dp)<1.0e-12_dp.and.abs(par(2)-1.0_dp)<1.0e-12_dp,'DPERKS family initial')
   if(fails/=0)error stop 1
   print *,'test_estimators: PASS'
contains
   subroutine chk(ok,name)
      logical,intent(in)::ok;character(*),intent(in)::name
      if(.not.ok)then;print *,trim(name),' failed';fails=fails+1;end if
   end subroutine
end program
