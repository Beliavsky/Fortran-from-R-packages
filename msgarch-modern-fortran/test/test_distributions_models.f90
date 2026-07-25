! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_distributions_models
   use msgarch
   use test_helpers
   implicit none
   character(len=8), parameter :: dists(6)=[character(len=8)::'norm','std','ged','snorm','sstd','sged']
   character(len=12), parameter :: models(5)=[character(len=12)::'sARCH','sGARCH','eGARCH','gjrGARCH','tGARCH']
   real(dp), parameter :: probs(3)=[0.05_dp,0.5_dp,0.95_dp]
   real(dp)::shape,skew,integral,h,z,w,q,p,meanv,varv,x,target,expected
   real(dp),allocatable::theta(:),lower(:),upper(:),draw(:),sorted_theta(:),unc(:)
   type(msgarch_spec)::spec,spec2,homogeneous_spec
   type(simulation_result)::simulation
   type(filter_result)::filtered
   logical::ok
   integer::i,j,n

   call seed_rng(12031)
   do i=1,size(dists)
      shape=8.0_dp;if(index(dists(i),'ged')>0)shape=1.6_dp
      skew=1.25_dp
      n=8000;h=80.0_dp/real(n,dp);integral=0.0_dp
      do j=0,n
         z=-40.0_dp+real(j,dp)*h
         if(j==0.or.j==n)then;w=1.0_dp;else if(mod(j,2)==0)then;w=2.0_dp;else;w=4.0_dp;end if
         integral=integral+w*innovation_pdf(z,dists(i),shape,skew)
      end do
      integral=integral*h/3.0_dp
      call assert_close(integral,1.0_dp,2.5e-3_dp,'density integrates to one: '//trim(dists(i)))
      do j=1,3
         p=probs(j)
         q=innovation_quantile(p,dists(i),shape,skew)
         call assert_close(innovation_cdf(q,dists(i),shape,skew),p,2.0e-7_dp,'cdf quantile inversion')
      end do
      allocate(draw(12000))
      do j=1,size(draw);draw(j)=random_innovation(dists(i),shape,skew);end do
      meanv=sum(draw)/real(size(draw),dp)
      varv=sum((draw-meanv)**2)/real(size(draw)-1,dp)
      call assert_true(abs(meanv)<0.08_dp,'innovation mean near zero: '//trim(dists(i)))
      call assert_true(abs(varv-1.0_dp)<0.14_dp,'innovation variance near one: '//trim(dists(i)))
      deallocate(draw)
   end do

   do i=1,size(models)
      spec=create_spec([models(i)],['sstd'])
      spec%regime(1)%shape=8.0_dp;spec%regime(1)%skew=1.1_dp
      call assert_true(spec_valid(spec),'valid single-regime model '//trim(models(i)))
      theta=pack_parameters(spec);call parameter_bounds(spec,lower,upper)
      call assert_true(size(theta)==parameter_count(spec),'parameter count')
      call assert_true(all(theta>=lower).and.all(theta<=upper),'default parameters in bounds')
      call unpack_parameters(spec,theta,spec2,ok);call assert_true(ok,'parameter round trip')
      call assert_close(maxval(abs(pack_parameters(spec2)-theta)),0.0_dp,1.0e-14_dp,'parameter pack round trip')
      simulation=simulate_msgarch(spec,180,2)
      call assert_all_finite(simulation%draw,'model simulation finite')
      filtered=hamilton_filter(spec,simulation%draw(1,:))
      call assert_true(filtered%loglik>-huge(1.0_dp)/2.0_dp,'model filter finite')
      call assert_true(all(filtered%variance>0.0_dp),'positive model variances')
      target=2.0_dp
      select case(trim(models(i)))
      case('sARCH');expected=target*(1.0_dp-spec%regime(1)%alpha)
      case('sGARCH');expected=target*(1.0_dp-spec%regime(1)%alpha-spec%regime(1)%beta)
      case('gjrGARCH');expected=target*(1.0_dp-spec%regime(1)%alpha-0.5_dp*spec%regime(1)%gamma-spec%regime(1)%beta)
      case('eGARCH');expected=log(target)*(1.0_dp-spec%regime(1)%beta)
      case('tGARCH');expected=target*(1.0_dp+0.5_dp*(spec%regime(1)%alpha+spec%regime(1)%gamma)-spec%regime(1)%beta)
      end select
      call assert_close(variance_target_intercept(target,spec%regime(1)),expected,1.0e-14_dp,'variance targeting formula')
   end do

   spec=create_spec([character(len=12)::'sGARCH','gjrGARCH'],[character(len=8)::'norm','std'],.true.)
   call assert_true(spec_valid(spec),'mixture specification valid')
   call assert_close(sum(stationary_distribution(spec%transition)),1.0_dp,1.0e-12_dp,'mixture stationary probabilities')
   x=bounded_map(0.3_dp,-2.0_dp,4.0_dp)
   call assert_close(bounded_unmap(x,-2.0_dp,4.0_dp),0.3_dp,1.0e-12_dp,'bounded map inversion')
   call assert_close(sum(simplex_mapping(simplex_unmapping([0.2_dp,0.3_dp,0.5_dp]))),1.0_dp,1.0e-12_dp,'simplex map')

   homogeneous_spec=create_spec([character(len=12)::'sGARCH','sGARCH'],[character(len=8)::'norm','norm'])
   homogeneous_spec%regime(1)%omega=0.20_dp;homogeneous_spec%regime(1)%alpha=0.10_dp;homogeneous_spec%regime(1)%beta=0.80_dp
   homogeneous_spec%regime(2)%omega=0.04_dp;homogeneous_spec%regime(2)%alpha=0.10_dp;homogeneous_spec%regime(2)%beta=0.80_dp
   homogeneous_spec%transition=reshape([0.90_dp,0.10_dp,0.20_dp,0.80_dp],[2,2],order=[2,1])
   sorted_theta=sort_parameters_by_variance(homogeneous_spec,pack_parameters(homogeneous_spec))
   call unpack_parameters(homogeneous_spec,sorted_theta,spec2,ok);unc=unconditional_variances(spec2)
   call assert_true(homogeneous_regimes(homogeneous_spec),'homogeneous regime detection')
   call assert_true(unc(1)<=unc(2),'regime parameters sorted by unconditional variance')
   call assert_close(spec2%transition(1,1),0.80_dp,1.0e-14_dp,'transition matrix permuted with regimes')

   write(*,'(a)')'Distribution and regime-model tests passed.'
end program test_distributions_models
