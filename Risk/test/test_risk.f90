! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
program test_risk
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use risk
   implicit none

   type(normal_distribution) :: normal
   type(uniform_distribution) :: uniform
   type(lognormal_distribution) :: lognormal
   type(exponential_distribution) :: exponential
   type(logistic_distribution) :: logistic
   type(student_t_distribution) :: student
   type(callback_distribution) :: triangular
   real(dp), parameter :: inf = huge(1.0_dp)
   real(dp), parameter :: z90 = 1.2815515655446004_dp
   real(dp), parameter :: es90 = -0.19499814659165204_dp
   real(dp), parameter :: mad_normal = 0.7978845608028654_dp
   real(dp) :: alpha(3), values(3)
   integer :: failures

   failures = 0
   normal = normal_distribution(0.0_dp,1.0_dp)
   uniform = uniform_distribution(0.0_dp,1.0_dp)
   lognormal = lognormal_distribution(0.0_dp,0.5_dp)
   exponential = exponential_distribution(2.0_dp)
   logistic = logistic_distribution(0.0_dp,1.0_dp)
   student = student_t_distribution(5.0_dp,0.0_dp,1.0_dp)
   call triangular%initialize(triangular_pdf,triangular_cdf,triangular_quantile,0.0_dp,1.0_dp)

   call assert_close('normal pdf',normal%pdf(0.0_dp),0.3989422804014327_dp,2.0e-13_dp,failures)
   call assert_close('normal cdf',normal%cdf(0.0_dp),0.5_dp,2.0e-15_dp,failures)
   call assert_close('normal quantile',normal%quantile(0.9_dp),z90,3.0e-12_dp,failures)
   call assert_close('lognormal inversion',lognormal%cdf(lognormal%quantile(0.8_dp)),0.8_dp,2.0e-12_dp,failures)
   call assert_close('uniform inversion',uniform%cdf(uniform%quantile(0.37_dp)),0.37_dp,2.0e-15_dp,failures)
   call assert_close('exponential inversion',exponential%cdf(exponential%quantile(0.72_dp)),0.72_dp,2.0e-15_dp,failures)
   call assert_close('logistic inversion',logistic%cdf(logistic%quantile(0.61_dp)),0.61_dp,2.0e-15_dp,failures)
   call assert_close('student t quantile',student%quantile(0.975_dp),2.57058183563631_dp,2.0e-9_dp,failures)
   call assert_close('student t inversion',student%cdf(student%quantile(0.91_dp)),0.91_dp,2.0e-11_dp,failures)
   call assert_close('callback quantile',varg(triangular,0.25_dp),0.5_dp,2.0e-15_dp,failures)
   call assert_close('callback expectation',expect(triangular,0.0_dp,1.0_dp),2.0_dp/3.0_dp,2.0e-10_dp,failures)

   call assert_close('varg',varg(normal,0.9_dp),z90,3.0e-12_dp,failures)
   call assert_close('esg',esg(normal,0.9_dp),es90,2.0e-8_dp,failures)
   call assert_close('tcm',tcm(normal,0.9_dp),1.6448536269514722_dp,4.0e-12_dp,failures)
   call assert_close('expp symmetry',expp(normal,0.5_dp,-inf,inf),0.0_dp,2.0e-7_dp,failures)
   call assert_close('expp uniform',expp(uniform,0.8_dp,0.0_dp,1.0_dp),0.38648820956430935_dp,2.0e-7_dp,failures)
   call assert_close('bvar',bvar(normal,0.9_dp,-inf),es90,2.0e-8_dp,failures)
   call assert_close('epsg',epsg(normal,0.9_dp),-0.01521578622618957_dp,2.0e-9_dp,failures)
   call assert_close('expect',expect(normal,-inf,inf),0.0_dp,2.0e-10_dp,failures)
   call assert_close('expvar',expvar(normal,0.9_dp,-inf,inf),0.9_dp,2.0e-8_dp,failures)
   call assert_close('omegag',omegag(normal,0.0_dp,-inf,inf),1.0_dp,2.0e-8_dp,failures)
   call assert_close('sortinog',sortinog(normal,0.0_dp,-inf,inf),0.0_dp,2.0e-10_dp,failures)
   call assert_close('kappag',kappag(normal,0.0_dp,3.0_dp,-inf,inf),0.0_dp,2.0e-10_dp,failures)
   call assert_close('wangg1',wangg1(uniform,0.5_dp,0.0_dp,1.0_dp),1.0_dp/6.0_dp,2.0e-9_dp,failures)
   call assert_close('wangg2',wangg2(uniform,0.5_dp,0.0_dp,1.0_dp),1.0_dp/6.0_dp,2.0e-9_dp,failures)
   call assert_close('stoneg1',stoneg1(uniform,0.0_dp,2.0_dp,0.0_dp,1.0_dp),1.0_dp/3.0_dp,2.0e-10_dp,failures)
   call assert_close('stoneg2',stoneg2(uniform,0.0_dp,2.0_dp,0.0_dp,1.0_dp),sqrt(1.0_dp/3.0_dp),2.0e-10_dp,failures)
   call assert_close('luceg1',luceg1(uniform,0.0_dp,1.0_dp,1.0_dp,0.0_dp),0.0_dp,2.0e-12_dp,failures)
   call assert_close('luceg2',luceg2(uniform,0.0_dp,1.0_dp,1.0_dp,0.0_dp),1.0_dp,2.0e-10_dp,failures)
   call assert_close('luceg3',luceg3(uniform,0.0_dp,1.0_dp,1.0_dp,0.0_dp),-1.0_dp,2.0e-7_dp,failures)
   call assert_close('luceg4',luceg4(uniform,0.0_dp,1.0_dp,1.0_dp,0.0_dp),1.0_dp,2.0e-10_dp,failures)
   call assert_close('saring1',saring1(normal,-inf,inf,1.0_dp,0.0_dp),1.0_dp,2.0e-9_dp,failures)
   call assert_close('saring2',saring2(uniform,0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp),1.0_dp,2.0e-10_dp,failures)
   call assert_close('saring3',saring3(uniform,0.0_dp,1.0_dp,2.0_dp,1.0_dp,0.0_dp),22.0_dp/45.0_dp,2.0e-9_dp,failures)
   call assert_close('bkg1',bkg1(normal,0.9_dp,-inf,inf),z90,2.0e-8_dp,failures)
   call assert_close('bkg2',bkg2(normal,0.9_dp,-inf,inf),es90,2.0e-8_dp,failures)
   call assert_close('bkg3',bkg3(normal,0.9_dp,-inf,inf,1.0_dp),es90-mad_normal,3.0e-8_dp,failures)
   call assert_close('bkg4',bkg4(normal,0.9_dp,-inf,inf,1.0_dp),z90+es90-mad_normal,3.0e-8_dp,failures)

   alpha = [0.1_dp,0.5_dp,0.9_dp]
   values = varg(normal,alpha)
   call assert_close('vector var first',values(1),-z90,4.0e-12_dp,failures)
   call assert_close('vector var middle',values(2),0.0_dp,2.0e-15_dp,failures)
   call assert_close('vector var last',values(3),z90,4.0e-12_dp,failures)
   values = expvar(normal,alpha,-inf,inf)
   call assert_close('vector expvar',values(2),0.5_dp,2.0e-8_dp,failures)
   values = esg(normal,alpha)
   call assert_close('vector esg',values(3),es90,2.0e-8_dp,failures)
   values = tcm(normal,alpha)
   call assert_close('vector tcm',values(3),1.6448536269514722_dp,4.0e-12_dp,failures)
   values = expp(uniform,[0.2_dp,0.5_dp,0.8_dp],0.0_dp,1.0_dp)
   call assert_close('vector expp',values(3),0.38648820956430935_dp,2.0e-7_dp,failures)
   values = bvar(normal,alpha,-inf)
   call assert_close('vector bvar',values(3),es90,2.0e-8_dp,failures)
   values = epsg(normal,[0.2_dp,0.8_dp,0.9_dp])
   call assert_close('vector epsg',values(3),-0.01521578622618957_dp,2.0e-9_dp,failures)
   values = omegag(normal,[0.0_dp,0.5_dp,1.0_dp],-inf,inf)
   call assert_close('vector omegag',values(1),1.0_dp,2.0e-8_dp,failures)
   values = sortinog(normal,[0.0_dp,0.5_dp,1.0_dp],-inf,inf)
   call assert_close('vector sortinog',values(1),0.0_dp,2.0e-10_dp,failures)
   values = kappag(normal,[0.0_dp,0.5_dp,1.0_dp],3.0_dp,-inf,inf)
   call assert_close('vector kappag',values(1),0.0_dp,2.0e-10_dp,failures)
   values = wangg1(uniform,[0.2_dp,0.5_dp,0.8_dp],0.0_dp,1.0_dp)
   call assert_close('vector wangg1',values(2),1.0_dp/6.0_dp,2.0e-9_dp,failures)
   values = wangg2(uniform,[0.2_dp,0.5_dp,0.8_dp],0.0_dp,1.0_dp)
   call assert_close('vector wangg2',values(2),1.0_dp/6.0_dp,2.0e-9_dp,failures)
   values = bkg1(normal,alpha,-inf,inf)
   call assert_close('vector bkg1',values(3),z90,2.0e-8_dp,failures)
   values = bkg2(normal,alpha,-inf,inf)
   call assert_close('vector bkg2',values(3),es90,2.0e-8_dp,failures)
   values = bkg3(normal,alpha,-inf,inf,1.0_dp)
   call assert_close('vector bkg3',values(3),es90-mad_normal,3.0e-8_dp,failures)
   values = bkg4(normal,alpha,-inf,inf,1.0_dp)
   call assert_close('vector bkg4',values(3),z90+es90-mad_normal,3.0e-8_dp,failures)

   call assert_nan('invalid alpha',varg(normal,1.1_dp),failures)
   call assert_nan('saring3 singular aa',saring3(uniform,0.0_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp),failures)

   if (failures > 0) then
      write(*,'(a,i0)') 'FAILED tests: ',failures
      error stop 1
   end if
   write(*,'(a)') 'All Risk tests passed.'

contains

   subroutine assert_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance
      integer, intent(inout) :: failures
      if (ieee_is_nan(actual)) then
         failures = failures+1
         write(*,'(a,1x,a)') trim(name),'returned NaN'
      else if (abs(actual-expected) > tolerance) then
         failures = failures+1
         write(*,'(a,2(1x,es24.16),1x,a,1x,es12.4)') trim(name),actual,expected,'tol',tolerance
      end if
   end subroutine assert_close

   subroutine assert_nan(name, actual, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual
      integer, intent(inout) :: failures
      if (.not. ieee_is_nan(actual)) then
         failures = failures+1
         write(*,'(a,1x,es24.16)') trim(name),actual
      end if
   end subroutine assert_nan

   function triangular_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x < 0.0_dp .or. x > 1.0_dp) then
         y = 0.0_dp
      else
         y = 2.0_dp*x
      end if
   end function triangular_pdf

   function triangular_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else if (x >= 1.0_dp) then
         p = 1.0_dp
      else
         p = x*x
      end if
   end function triangular_cdf

   function triangular_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x
      x = sqrt(max(0.0_dp,min(1.0_dp,p)))
   end function triangular_quantile

end program test_risk
