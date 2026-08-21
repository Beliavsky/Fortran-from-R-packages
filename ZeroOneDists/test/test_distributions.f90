program test_distributions
   use zero_one_dists
   implicit none
   integer :: fails
   real(dp) :: base_mu, theta
   fails = 0
   call check(dber(0.37_dp,0.42_dp,3.7_dp,0.18_dp),1.3994300718592987_dp,2.0e-12_dp,'dBER',fails)
   call check(pber(0.37_dp,0.42_dp,3.7_dp,0.18_dp),0.4351276620721301_dp,2.0e-12_dp,'pBER',fails)
   call check(qber(0.73_dp,0.42_dp,3.7_dp,0.18_dp),0.5982181243621164_dp,2.0e-11_dp,'qBER',fails)
   call ber2_transform(0.63_dp,0.54_dp,base_mu,theta)
   call check(base_mu,0.7165223184543638_dp,2.0e-14_dp,'BER2 base mu',fails)
   call check(theta,0.3996_dp,2.0e-14_dp,'BER2 theta',fails)
   call check(dber2(0.37_dp,0.63_dp,2.8_dp,0.54_dp),0.7462487720893466_dp,2.0e-12_dp,'dBER2',fails)
   call check(pber2(0.37_dp,0.63_dp,2.8_dp,0.54_dp),0.20957088955819678_dp,2.0e-12_dp,'pBER2',fails)
   call check(qber2(0.73_dp,0.63_dp,2.8_dp,0.54_dp),0.8590531993395591_dp,2.0e-11_dp,'qBER2',fails)
   call check(duhlg(0.37_dp,1.4_dp),1.0642810553410946_dp,2.0e-13_dp,'dUHLG',fails)
   call check(puhlg(0.37_dp,1.4_dp),0.4562268803945746_dp,2.0e-13_dp,'pUHLG',fails)
   call check(quhlg(0.73_dp,1.4_dp),0.6542893725992317_dp,2.0e-13_dp,'qUHLG',fails)
   call check(dumb(0.37_dp,0.8_dp),1.9233490647660008_dp,3.0e-13_dp,'dUMB',fails)
   call check(pumb(0.37_dp,0.8_dp),0.6720178763726078_dp,3.0e-13_dp,'pUMB',fails)
   call check(qumb(0.73_dp,0.8_dp),0.40217863796399383_dp,2.0e-12_dp,'qUMB',fails)
   call check(duphn(0.37_dp,1.7_dp,0.9_dp),1.7627477957142976_dp,5.0e-12_dp,'dUPHN',fails)
   call check(puphn(0.37_dp,1.7_dp,0.9_dp),0.09741039484872549_dp,5.0e-13_dp,'pUPHN',fails)
   call check(quphn(0.73_dp,1.7_dp,0.9_dp),0.6428742066930911_dp,2.0e-11_dp,'qUPHN',fails)
   if (fails /= 0) error stop 1
   print '(a)', 'test_distributions: PASS'
contains
   subroutine check(got,want,tol,name,nfail)
      real(dp), intent(in) :: got,want,tol
      character(*), intent(in) :: name
      integer, intent(inout) :: nfail
      if (abs(got-want) > tol*max(1.0_dp,abs(want))) then
         print '(a,2(1x,es24.16))', trim(name)//' FAIL:',got,want
         nfail=nfail+1
      end if
   end subroutine check
end program test_distributions
