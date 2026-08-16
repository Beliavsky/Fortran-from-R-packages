program test_family
   use discretedists
   implicit none
   type(discrete_family_t)::f
   real(dp)::s1,s2,c11,c12,c22,m,v,eta,x,r
   integer::fails
   fails=0
   f=berg();call f%score(2.0_dp,2.0_dp,2.2_dp,s1,s2)
   call f%curvature(2.0_dp,2.0_dp,2.2_dp,c11,c12,c22)
   call chk(abs(s1)<100.0_dp.and.abs(s2)<100.0_dp,'BerG score')
   call chk(abs(c11)<1.0e6_dp.and.abs(c12)<1.0e6_dp.and.abs(c22)<1.0e6_dp,'BerG curvature')
   call f%mean_variance(2.0_dp,2.2_dp,m,v)
   call chk(abs(m-2.0_dp)<1.0e-12_dp.and.abs(v-4.4_dp)<1.0e-12_dp,'BerG moments')
   f=ggeo();eta=f%linkfun(0.4_dp);x=f%linkinv(eta)
   call chk(abs(x-0.4_dp)<1.0e-12_dp,'GGEO link')
   call f%score(2.0_dp,0.4_dp,1.3_dp,s1,s2)
   call chk(abs(s1)<100.0_dp.and.abs(s2)<100.0_dp,'GGEO score')
   r=f%rqres(2,0.4_dp,1.3_dp,0.5_dp)
   call chk(abs(r)<10.0_dp,'rqres')
   f=compo();call f%mean_variance(2.0_dp,1.0_dp,m,v)
   call chk(abs(m-2.0_dp)<1.0e-8_dp.and.abs(v-2.0_dp)<1.0e-8_dp,'COMPO moments')
   call chk(f%valid(2.0_dp,1.0_dp),'COMPO valid')
   f=ggeo()
   call chk(.not.f%valid(1.2_dp,1.0_dp),'GGEO invalid')
   if(fails/=0)error stop 1
   print *,'test_family: PASS'
contains
   subroutine chk(ok,name)
      logical,intent(in)::ok;character(*),intent(in)::name
      if(.not.ok)then;print *,trim(name),' failed';fails=fails+1;end if
   end subroutine
end program
