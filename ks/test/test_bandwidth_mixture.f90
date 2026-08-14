! SPDX-License-Identifier: GPL-2.0-only
program test_bandwidth_mixture
  use ks, only: dp, hns_1d, hns_matrix, hns_diag, hlscv_1d, hpi, hscv, hpi_matrix, &
                mise_normal_mixture, amise_normal_mixture, normal_mixture_modes, determinant_spd
  implicit none
  real(dp)::x1(9),x2(9,2),bw,bw2,bw3,Hm(2,2),Hd(2,2),Hp(2,2),mus(2,2),sig(2,2,2),props(2), &
            mise,amise,modes(2,2),det_hm,det_hd
  integer::i
  x1=[-1.4_dp,-0.9_dp,-0.4_dp,-0.1_dp,0.2_dp,0.5_dp,0.9_dp,1.5_dp,2.0_dp]
  do i=1,9;x2(i,1)=x1(i);x2(i,2)=0.4_dp*x1(i)+0.25_dp*sin(real(i,dp));end do
  bw=hns_1d(x1);bw2=hpi(x1);bw3=hscv(x1)
  if(min(bw,bw2,bw3)<=0.0_dp) error stop '1d bandwidth'
  if(hlscv_1d(x1)<=0.0_dp) error stop 'lscv bandwidth'
  call hns_matrix(x2,Hm);call hns_diag(x2,Hd)
  det_hm=determinant_spd(Hm);det_hd=determinant_spd(Hd)
  if(det_hm<=0.0_dp.or.det_hd<=0.0_dp) error stop 'Hns SPD'
  call hpi_matrix(x2,Hp,maxiter=250)
  if(determinant_spd(Hp)<=0.0_dp) error stop 'Hpi SPD'
  mus(1,:)=[-1.0_dp,0.0_dp];mus(2,:)=[1.2_dp,0.8_dp];props=[0.4_dp,0.6_dp]
  sig=0.0_dp;sig(:,:,1)=reshape([0.5_dp,0.1_dp,0.1_dp,0.8_dp],[2,2]);sig(:,:,2)=reshape([0.9_dp,-0.15_dp,-0.15_dp,0.6_dp],[2,2])
  mise=mise_normal_mixture(Hm,mus,sig,props,100);amise=amise_normal_mixture(Hm,mus,sig,props,100)
  if(.not.(mise>=0.0_dp).or..not.(amise>=0.0_dp)) error stop 'mixture risk'
  call normal_mixture_modes(mus,sig,props,modes,maxiter=300)
  if(any(abs(modes)>10.0_dp)) error stop 'mixture modes'
  print *, 'test_bandwidth_mixture: PASS'
end program
