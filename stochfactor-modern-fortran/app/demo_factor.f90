! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program demo_factor
  use sv_kinds,only:dp
  use sv_rng,only:seed_rng
  use fsv_types
  use fsv_core
  implicit none
  real(dp)::b(4,2),ip(4,3),fp(2,2),cov(4,4)
  type(fsv_sim_result)::sim
  type(fsv_options)::opt
  type(fsv_draws)::draws
  integer::i
  call seed_rng(2027);b=0.0_dp;b(:,1)=[1.0_dp,.8_dp,.5_dp,.2_dp];b(2:,2)=[1.0_dp,.4_dp,-.3_dp]
  do i=1,4;ip(i,:)=[-1.5_dp+.1_dp*i,.9_dp+.01_dp*i,.2_dp];end do;fp=reshape([.97_dp,.95_dp,.12_dp,.18_dp],[2,2])
  call simulate_fsv(150,b,ip,fp,sim);opt%draws=40;opt%burnin=20;call fit_fsv(sim%y,2,opt,draws);call running_covariance(draws,150,cov)
  print '(a)','Posterior mean final covariance matrix:'
  do i=1,4;print '(4f11.5)',cov(i,:);end do
end program demo_factor
