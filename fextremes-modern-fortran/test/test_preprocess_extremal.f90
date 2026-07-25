! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software under GPL-2.0-or-later.
program test_preprocess_extremal
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state,seed_rng
  use fextremes_preprocess
  use fextremes_extremal_index
  implicit none
  real(dp)::x(10),sim(12000),probs(2)
  real(dp),allocatable::mx(:)
  integer,allocatable::mi(:)
  type(point_process_result)::pp
  type(decluster_result)::dc
  type(theta_result)::tr
  type(rng_state)::rng
  x=[1.0_dp,4.0_dp,2.0_dp,3.0_dp,9.0_dp,5.0_dp,7.0_dp,6.0_dp,8.0_dp,0.0_dp]
  call block_maxima(x,3,mx,mi)
  call assert_true(all(abs(mx-[4.0_dp,9.0_dp,8.0_dp,0.0_dp])<1.0e-12_dp),'block maxima')
  call assert_true(all(mi==[2,5,9,10]),'block indices')
  call assert_true(abs(find_threshold(x,3)-6.0_dp)<1.0e-12_dp,'find threshold')
  call point_process(x,6.0_dp,pp)
  call assert_true(all(pp%indices==[5,7,9]),'point process indices')
  call decluster(pp%values,pp%indices,2,dc)
  call assert_true(size(dc%maxima)==1 .and. abs(dc%maxima(1)-9.0_dp)<1.0e-12_dp,'declustering')
  call seed_rng(rng,1985); call theta_simulate(rng,'max',size(sim),0.5_dp,sim); probs=[0.95_dp,0.98_dp]
  call cluster_theta(sim,22,probs,tr)
  call assert_true(all(tr%theta>0.25_dp .and. tr%theta<0.8_dp),'cluster theta')
  call ferro_seg_theta(sim,probs,tr)
  call assert_true(all(tr%theta>0.2_dp .and. tr%theta<=1.0_dp),'Ferro-Segers theta')
  print '(a)','Preprocessing and extremal-index tests passed.'
contains
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond; character(len=*),intent(in)::msg
    if(.not.cond) then; print *,trim(msg); error stop 1; end if
  end subroutine
end program test_preprocess_extremal
