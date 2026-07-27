! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_univariate_extended
   use gogarch
   use test_helpers
   implicit none
   integer, parameter :: n = 700, nsim = 60000
   character(len=5), parameter :: distributions(6) = ['norm ','snorm','std  ','sstd ','ged  ','sged ']
   real(dp) :: y(n), h(n), power(n), residuals(n), filtered_power(n), variance(n), standardized(n)
   real(dp) :: mf(5), vf(5), samples(nsim), mean_sample, variance_sample, density_sum, step, x
   real(dp) :: shape, skew, ll
   real(dp), allocatable :: arch(:), leverage(:), garch(:)
   type(univariate_spec) :: spec
   type(garch11_fit) :: fit
   logical :: ok
   integer :: i, j

   call seed_rng(112233)
   do j = 1, size(distributions)
      select case (trim(distributions(j)))
      case ('std','sstd')
         shape = 7.5_dp
      case ('ged','sged')
         shape = 1.35_dp
      case default
         shape = 8.0_dp
      end select
      if (distributions(j)(1:1) == 's') then
         skew = 1.45_dp
      else
         skew = 1.0_dp
      end if
      do i = 1, nsim
         samples(i) = random_innovation(trim(distributions(j)),shape,skew)
      end do
      mean_sample = sum(samples)/real(nsim,dp)
      variance_sample = sum((samples-mean_sample)**2)/real(nsim,dp)
      call assert_true(abs(mean_sample) < 0.04_dp,trim(distributions(j))//' RNG mean')
      call assert_true(abs(variance_sample-1.0_dp) < 0.07_dp,trim(distributions(j))//' RNG variance')
      step = 0.01_dp
      density_sum = 0.0_dp
      do i = -4000, 4000
         x = step*real(i,dp)
         density_sum = density_sum+innovation_pdf(x,trim(distributions(j)),shape,skew)
      end do
      density_sum = density_sum*step
      call assert_true(abs(density_sum-1.0_dp) < 0.015_dp,trim(distributions(j))//' density normalization')
      call simulate_garchpq(350,0.0_dp,0.04_dp,[0.08_dp],[0.86_dp],trim(distributions(j)),shape,skew, &
         y(1:350),h(1:350),burnin=400)
      spec = univariate_spec()
      spec%distribution = trim(distributions(j))
      spec%shape = shape
      spec%skew = skew
      spec%fit_shape = .false.
      spec%fit_skew = .false.
      fit = fit_univariate(y(1:350),spec,max_iterations=320)
      call assert_true(fit%status <= 1,trim(distributions(j))//' likelihood fit')
      call assert_true(trim(fit%distribution) == trim(distributions(j)),trim(distributions(j))//' fit label')
   end do

   allocate(arch(2),garch(2),leverage(0))
   arch = [0.06_dp,0.03_dp]
   garch = [0.70_dp,0.15_dp]
   call simulate_garchpq(n,0.01_dp,0.03_dp,arch,garch,'std',8.0_dp,1.0_dp,y,h,burnin=600)
   call filter_garchpq(y,0.01_dp,0.03_dp,arch,garch,'std',8.0_dp,1.0_dp,residuals,variance,standardized,ll,ok)
   call assert_true(ok,'GARCH(2,2) Student filter')
   call assert_all_finite(standardized,'GARCH(2,2) standardized residuals')
   spec%model = 'garch'
   spec%distribution = 'std'
   spec%p = 2
   spec%q = 2
   spec%shape = 8.0_dp
   spec%fit_shape = .true.
   fit = fit_univariate(y,spec,max_iterations=650)
   call assert_true(fit%status <= 1,'GARCH(2,2) Student fit')
   call assert_true(size(fit%arch) == 2 .and. size(fit%garch) == 2,'higher-order coefficient retention')
   call assert_true(fit%shape > 2.0_dp,'fitted Student shape')
   call forecast_univariate(fit,5,mf,vf)
   call assert_true(all(vf > 0.0_dp),'GARCH(2,2) forecasts')

   deallocate(arch,garch,leverage)
   allocate(arch(1),leverage(1),garch(1))
   arch = [0.08_dp]
   leverage = [0.22_dp]
   garch = [0.84_dp]
   call simulate_aparch(n,0.0_dp,0.035_dp,arch,leverage,garch,1.35_dp,'sged',1.45_dp,1.30_dp,y,h, &
      burnin=700,power_scale=power)
   call filter_aparch(y,0.0_dp,0.035_dp,arch,leverage,garch,1.35_dp,'sged',1.45_dp,1.30_dp,residuals, &
      filtered_power,variance,standardized,ll,ok)
   call assert_true(ok,'APARCH skew-GED filter')
   call assert_close(filtered_power(2),0.035_dp+0.08_dp*(abs(residuals(1))-0.22_dp*residuals(1))**1.35_dp+ &
      0.84_dp*filtered_power(1),1.0e-12_dp,'APARCH filter recursion')
   spec%model = 'aparch'
   spec%distribution = 'sged'
   spec%p = 1
   spec%o = 1
   spec%q = 1
   spec%delta = 1.5_dp
   spec%shape = 1.5_dp
   spec%skew = 1.2_dp
   spec%fit_delta = .true.
   spec%fit_shape = .true.
   spec%fit_skew = .true.
   fit = fit_univariate(y,spec,max_iterations=850)
   call assert_true(fit%status <= 1,'APARCH skew-GED fit')
   call assert_true(fit%delta > 0.25_dp .and. fit%delta < 4.0_dp,'fitted APARCH delta')
   call assert_true(abs(fit%leverage(1)) < 1.0_dp,'fitted APARCH leverage')
   call assert_true(fit%shape > 0.25_dp .and. fit%skew > 0.2_dp,'fitted skew-GED parameters')
   call forecast_univariate(fit,5,mf,vf)
   call assert_true(all(vf > 0.0_dp),'APARCH forecasts')

   write(*,'(a)') 'Higher-order GARCH, APARCH, and distribution tests passed.'
end program test_univariate_extended
