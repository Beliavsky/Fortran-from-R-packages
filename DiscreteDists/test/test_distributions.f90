program test_distributions
   use discretedists
   implicit none
   real(dp),parameter :: tol=5.0e-7_dp
   integer :: fails
   fails=0

   call close(dcompo(3.0_dp,2.0_dp,1.0_dp),poisson3(2.0_dp),tol,'COMPO Poisson pmf')
   call close(pcompo(3.0_dp,2.0_dp,1.0_dp),poisson_cdf3(2.0_dp),tol,'COMPO Poisson cdf')
   call close(dcompo2(3.0_dp,2.0_dp,0.0_dp),dcompo(3.0_dp,2.0_dp,1.0_dp),tol,'COMPO2 phi=0')
   call close(dhyperpo(3.0_dp,2.0_dp,1.0_dp),dcompo(3.0_dp,2.0_dp,1.0_dp),tol,'HYPERPO sigma=1')
   call close(dhyperpo2(3.0_dp,2.0_dp,1.0_dp),dcompo(3.0_dp,2.0_dp,1.0_dp),tol,'HYPERPO2 sigma=1')
   call close(ddgeii(3.0_dp,0.4_dp,1.0_dp),0.6_dp*0.4_dp**3,tol,'DGEII geometric')
   call close(dggeo(3.0_dp,0.4_dp,1.0_dp),0.6_dp*0.4_dp**3,tol,'GGEO geometric')

   call cdf_sum_checks()
   call quantile_checks()

   if(fails/=0) error stop 1
   print *,'test_distributions: PASS'
contains
   real(dp) function poisson3(mu)
      real(dp),intent(in) :: mu
      poisson3=exp(-mu)*mu**3/6.0_dp
   end function poisson3

   real(dp) function poisson_cdf3(mu)
      real(dp),intent(in) :: mu
      poisson_cdf3=exp(-mu)*(1.0_dp+mu+mu**2/2.0_dp+mu**3/6.0_dp)
   end function poisson_cdf3

   subroutine close(a,b,t,name)
      real(dp),intent(in) :: a,b,t
      character(*),intent(in) :: name
      if(abs(a-b)>t*max(1.0_dp,abs(a),abs(b))) then
         print *,trim(name),' failed:',a,b
         fails=fails+1
      end if
   end subroutine close

   subroutine cdf_sum_checks()
      integer :: j
      real(dp) :: s
      s=0.0_dp
      do j=0,4; s=s+dberg(real(j,dp),2.0_dp,2.2_dp); end do
      call close(pberg(4.0_dp,2.0_dp,2.2_dp),s,tol,'BerG cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddbh(real(j,dp),0.45_dp); end do
      call close(pdbh(4.0_dp,0.45_dp),s,tol,'DBH cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddikum(real(j,dp),1.3_dp,2.1_dp); end do
      call close(pdikum(4.0_dp,1.3_dp,2.1_dp),s,tol,'DIKUM cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddld(real(j,dp),0.7_dp); end do
      call close(pdld(4.0_dp,0.7_dp),s,tol,'DLD corrected cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddmolbe(real(j,dp),1.4_dp,0.8_dp); end do
      call close(pdmolbe(4.0_dp,1.4_dp,0.8_dp),s,tol,'DMOLBE cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddperks(real(j,dp),0.8_dp,0.7_dp); end do
      call close(pdperks(4.0_dp,0.8_dp,0.7_dp),s,tol,'DPERKS cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+ddspa(real(j,dp),1.2_dp,0.6_dp); end do
      call close(pdspa(4.0_dp,1.2_dp,0.6_dp),s,tol,'DsPA cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+dggeo(real(j,dp),0.4_dp,1.3_dp); end do
      call close(pggeo(4.0_dp,0.4_dp,1.3_dp),s,tol,'GGEO cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+dhyperpo(real(j,dp),1.8_dp,1.4_dp); end do
      call close(phyperpo(4.0_dp,1.8_dp,1.4_dp),s,1.0e-9_dp,'HYPERPO cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+dhyperpo2(real(j,dp),1.8_dp,1.4_dp); end do
      call close(phyperpo2(4.0_dp,1.8_dp,1.4_dp),s,1.0e-9_dp,'HYPERPO2 cdf/pmf')
      s=0.0_dp
      do j=0,4; s=s+dpoisxl(real(j,dp),0.7_dp); end do
      call close(ppoisxl(4.0_dp,0.7_dp),s,tol,'POISXL cdf/pmf')
   end subroutine cdf_sum_checks

   subroutine qcheck(q,cdf,p,name)
      integer,intent(in) :: q
      real(dp),intent(in) :: cdf,p
      character(*),intent(in) :: name
      if(q<0 .or. cdf<p) then
         print *,trim(name),' quantile failed'
         fails=fails+1
      end if
   end subroutine qcheck

   subroutine quantile_checks()
      real(dp),parameter :: p=0.73_dp
      integer :: q
      q=qberg(p,2.0_dp,2.2_dp)
      call qcheck(q,pberg(real(q,dp),2.0_dp,2.2_dp),p,'BerG')
      q=qcompo(p,2.0_dp,1.2_dp)
      call qcheck(q,pcompo(real(q,dp),2.0_dp,1.2_dp),p,'COMPO')
      q=qcompo2(p,2.0_dp,0.2_dp)
      call qcheck(q,pcompo2(real(q,dp),2.0_dp,0.2_dp),p,'COMPO2')
      q=qdbh(p,0.45_dp); call qcheck(q,pdbh(real(q,dp),0.45_dp),p,'DBH')
      q=qdgeii(p,0.4_dp,1.4_dp)
      call qcheck(q,pdgeii(real(q,dp),0.4_dp,1.4_dp),p,'DGEII')
      q=qdikum(p,1.3_dp,2.1_dp)
      call qcheck(q,pdikum(real(q,dp),1.3_dp,2.1_dp),p,'DIKUM')
      q=qdld(p,0.7_dp); call qcheck(q,pdld(real(q,dp),0.7_dp),p,'DLD')
      q=qdmolbe(p,1.4_dp,0.8_dp)
      call qcheck(q,pdmolbe(real(q,dp),1.4_dp,0.8_dp),p,'DMOLBE')
      q=qdperks(p,0.8_dp,0.7_dp)
      call qcheck(q,pdperks(real(q,dp),0.8_dp,0.7_dp),p,'DPERKS')
      q=qdspa(p,1.2_dp,0.6_dp)
      call qcheck(q,pdspa(real(q,dp),1.2_dp,0.6_dp),p,'DsPA')
      q=qggeo(p,0.4_dp,1.3_dp)
      call qcheck(q,pggeo(real(q,dp),0.4_dp,1.3_dp),p,'GGEO')
      q=qhyperpo(p,1.8_dp,1.4_dp)
      call qcheck(q,phyperpo(real(q,dp),1.8_dp,1.4_dp),p,'HYPERPO')
      q=qhyperpo2(p,1.8_dp,1.4_dp)
      call qcheck(q,phyperpo2(real(q,dp),1.8_dp,1.4_dp),p,'HYPERPO2')
      q=qpoisxl(p,0.7_dp)
      call qcheck(q,ppoisxl(real(q,dp),0.7_dp),p,'POISXL')
   end subroutine quantile_checks
end program test_distributions
