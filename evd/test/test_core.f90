program test_core
   use evd
   implicit none
   integer :: fails, i
   real(dp) :: p, q, back, mar1(3), mar2(3), depa(1), asy(1,3), mar3(3,3)
   real(dp) :: x3(3), probs(5), dat(12,2), epdat(12,2), a(3), xx(3), scale(2), shape(2), thr(2), nll
   type(cluster_result_t) :: cl
   fails=0

   probs=[0.05_dp,0.2_dp,0.5_dp,0.8_dp,0.95_dp]
   do i=1,size(probs)
      p=probs(i)
      q=qgev(p,1.2_dp,2.3_dp,0.2_dp)
      back=pgev(q,1.2_dp,2.3_dp,0.2_dp)
      call check(abs(back-p)<2e-12_dp,'GEV p/q positive shape',fails)
      q=qgev(p,-0.4_dp,1.1_dp,-0.25_dp)
      back=pgev(q,-0.4_dp,1.1_dp,-0.25_dp)
      call check(abs(back-p)<2e-12_dp,'GEV p/q negative shape',fails)
      q=qgpd(p,0.7_dp,1.8_dp,0.15_dp)
      back=pgpd(q,0.7_dp,1.8_dp,0.15_dp)
      call check(abs(back-p)<2e-12_dp,'GPD p/q',fails)
   end do

   mar1=[0.0_dp,1.0_dp,0.0_dp]
   mar2=[0.0_dp,1.0_dp,0.0_dp]
   p=pbvlog(0.4_dp,-0.2_dp,1.0_dp,mar1,mar2)
   back=pgev(0.4_dp,0.0_dp,1.0_dp,0.0_dp)*pgev(-0.2_dp,0.0_dp,1.0_dp,0.0_dp)
   call check(abs(p-back)<2e-12_dp,'bivariate logistic independence CDF',fails)
   call check(abs(abvlog(0.37_dp,1.0_dp)-1.0_dp)<2e-12_dp,'Pickands independence',fails)
   call check(dbvlog(0.4_dp,-0.2_dp,0.7_dp,mar1,mar2)>0.0_dp,'bivariate density positive',fails)

   mar3=0.0_dp
   mar3(:,2)=1.0_dp
   x3=[0.2_dp,-0.1_dp,0.7_dp]
   call check(pmvlog(x3,1.0_dp,mar3)>0.0_dp .and. pmvlog(x3,1.0_dp,mar3)<1.0_dp,'multivariate logistic CDF',fails)
   call check(dmvlog(x3,0.7_dp,mar3)>0.0_dp,'multivariate logistic density',fails)
   depa=[0.7_dp]
   asy=1.0_dp
   call check(abs(pmvalog(x3,depa,asy,mar3)-pmvlog(x3,0.7_dp,mar3))<1e-12_dp,'alog reduces to log',fails)

   dat(:,1)=[0.1_dp,0.4_dp,0.2_dp,0.9_dp,1.1_dp,0.5_dp,1.7_dp,0.3_dp,1.5_dp,0.7_dp,2.0_dp,1.3_dp]
   dat(:,2)=[0.2_dp,0.1_dp,0.6_dp,0.8_dp,0.4_dp,1.2_dp,1.5_dp,0.5_dp,1.7_dp,0.9_dp,1.8_dp,1.4_dp]
   call empirical_to_exponential(dat,epdat)
   xx=[0.0_dp,0.5_dp,1.0_dp]
   a=abvnonpar(xx,epdat,'cfg')
   call check(all(a>=max(xx,1.0_dp-xx)) .and. all(a<=1.0_dp),'nonparametric Pickands bounds',fails)

   thr=[0.5_dp,0.5_dp]
   scale=[1.0_dp,1.2_dp]
   shape=[0.1_dp,-0.05_dp]
   nll=bvpot_censored_nll(dat,thr,'log',[0.7_dp],scale,shape)
   call check(nll<1.0e90_dp,'censored POT likelihood finite',fails)
   nll=bvpot_poisson_nll(dat,thr,'log',[0.7_dp],scale,shape)
   call check(nll<1.0e90_dp,'Poisson POT likelihood finite',fails)

   cl=clusters([0.0_dp,2.0_dp,1.5_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp], &
               [1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],2, &
               [-huge(1.0_dp),-huge(1.0_dp),-huge(1.0_dp),-huge(1.0_dp),-huge(1.0_dp),-huge(1.0_dp),-huge(1.0_dp)],1)
   call check(cl%nclusters>=1,'clusters smoke',fails)

   if(fails>0) then
      write(*,'(a,i0)') 'test_core: FAIL ',fails
      error stop 1
   end if
   write(*,'(a)') 'test_core: PASS'
contains
   subroutine check(ok,msg,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      integer,intent(inout)::fails
      if(.not.ok) then
      write(*,'(a)') 'FAIL: '//msg
      fails=fails+1
      end if
   end subroutine
end program test_core
