! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
program fit_csv
  use sv_kinds,only:dp
  use sv_rng,only:seed_rng
  use sv_types
  use sv_core
  use fsv_types
  use fsv_core
  implicit none
  character(256)::file,mode,arg,line
  real(dp),allocatable::y(:,:),v(:)
  integer::n,m,ios,i,j,r
  call get_command_argument(1,file);call get_command_argument(2,mode);if(len_trim(file)==0)error stop 'usage: fit_csv file.csv sv|factor [factors]'
  open(10,file=trim(file),status='old',action='read');read(10,'(a)',iostat=ios)line;n=0
  do;read(10,'(a)',iostat=ios)line;if(ios/=0)exit;n=n+1;end do;rewind(10);read(10,'(a)')line
  m=1;do i=1,len_trim(line);if(line(i:i)==',')m=m+1;end do;m=m-1;allocate(y(n,m))
  do i=1,n;read(10,*,iostat=ios)arg,(y(i,j),j=1,m);if(ios/=0)error stop 'CSV parse error';end do;close(10);call seed_rng(13579)
  if(trim(mode)=='sv')then
    block
      type(sv_params)::p;type(sv_prior)::pr;type(sv_mcmc_options)::op;type(sv_draws)::dr
      allocate(v(n));v=y(:,1);p%mu=log(max(sum(v*v)/n,1.0e-6_dp));op%draws=100;op%burnin=50;call svsample(v,p,pr,op,dr);print '(a,3f12.6)','mu phi sigma ',sum(dr%para(1,:))/dr%ndraws,sum(dr%para(2,:))/dr%ndraws,sum(dr%para(3,:))/dr%ndraws
    end block
  else
    call get_command_argument(3,arg);read(arg,*,iostat=ios)r;if(ios/=0)r=1
    block
      type(fsv_options)::op;type(fsv_draws)::dr;real(dp),allocatable::cov(:,:)
      op%draws=50;op%burnin=25;call fit_fsv(y,r,op,dr);allocate(cov(m,m));call running_covariance(dr,n,cov);print '(a)','posterior mean final covariance';do i=1,m;print '(*(f12.6,1x))',cov(i,:);end do
    end block
  end if
end program fit_csv
