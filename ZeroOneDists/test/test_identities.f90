program test_identities
   use zero_one_dists
   implicit none
   integer :: fails, i
   real(dp), parameter :: probs(7)=[0.01_dp,0.1_dp,0.25_dp,0.5_dp,0.75_dp,0.9_dp,0.99_dp]
   real(dp) :: x, h, dnum, mean2
   fails=0
   do i=1,size(probs)
      x=qber(probs(i),0.35_dp,5.0_dp,0.2_dp)
      call close(pber(x,0.35_dp,5.0_dp,0.2_dp),probs(i),2.0e-11_dp,'BER inversion',fails)
      x=qber2(probs(i),0.61_dp,4.0_dp,0.5_dp)
      call close(pber2(x,0.61_dp,4.0_dp,0.5_dp),probs(i),2.0e-11_dp,'BER2 inversion',fails)
      x=quhlg(probs(i),0.7_dp)
      call close(puhlg(x,0.7_dp),probs(i),2.0e-13_dp,'UHLG inversion',fails)
      x=qumb(probs(i),1.2_dp)
      call close(pumb(x,1.2_dp),probs(i),3.0e-12_dp,'UMB inversion',fails)
      x=quphn(probs(i),1.3_dp,0.6_dp)
      call close(puphn(x,1.3_dp,0.6_dp),probs(i),3.0e-11_dp,'UPHN inversion',fails)
   end do
   h=1.0e-6_dp
   x=0.43_dp
   dnum=(puphn(x+h,1.3_dp,0.6_dp)-puphn(x-h,1.3_dp,0.6_dp))/(2.0_dp*h)
   call close(dnum,duphn(x,1.3_dp,0.6_dp),5.0e-8_dp,'UPHN d/dx CDF',fails)
   mean2=integral_mean_ber2(0.61_dp,4.0_dp,0.5_dp)
   call close(mean2,0.61_dp,2.0e-5_dp,'BER2 mean parameterization',fails)
   call close(pber2(0.37_dp,0.5_dp,2.0_dp,1.0_dp),0.37_dp,1.0e-14_dp,'BER2 uniform boundary',fails)
   call close(dber2(0.37_dp,0.5_dp,2.0_dp,1.0_dp),1.0_dp,1.0e-14_dp,'BER2 uniform density',fails)
   if(fails/=0)error stop 1
   print '(a)','test_identities: PASS'
contains
   real(dp) function integral_mean_ber2(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      integer,parameter::n=20000
      integer::j
      real(dp)::xx,hh,s
      hh=1.0_dp/real(n,dp)
      s=0.0_dp
      do j=1,n-1
         xx=real(j,dp)*hh
         if(mod(j,2)==0)then
            s=s+2.0_dp*xx*dber2(xx,mu,sigma,nu)
         else
            s=s+4.0_dp*xx*dber2(xx,mu,sigma,nu)
         end if
      end do
      v=hh*s/3.0_dp
   end function integral_mean_ber2
   subroutine close(got,want,tol,name,nfail)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::name
      integer,intent(inout)::nfail
      if(abs(got-want)>tol*max(1.0_dp,abs(want)))then
         print '(a,2(1x,es24.16))',trim(name)//' FAIL:',got,want
         nfail=nfail+1
      end if
   end subroutine close
end program test_identities
