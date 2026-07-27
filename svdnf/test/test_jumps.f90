! SPDX-License-Identifier: GPL-3.0-only
program test_jumps
  use svdnf
  use test_support
  implicit none
  type(svm_dynamics) :: bates, pmd, dps
  type(simulation_result) :: simulated
  type(filter_result) :: filtered
  type(grid_type) :: grids
  real(dp), allocatable :: transition(:,:)

  bates=dynamics_svm('Bates',mu=0.03_dp,kappa=2.5_dp,theta=0.04_dp,sigma=0.35_dp, &
    rho=-0.5_dp,omega=40.0_dp,alpha=-0.02_dp,delta=0.04_dp,h=1.0_dp/252.0_dp)
  simulated=model_simulate(bates,200,initial_volatility=0.04_dp,seed=999)
  call assert_true(simulated%ok,'Bates simulation')
  call assert_true(sum(simulated%jump_counts)>0,'Bates simulated jumps')
  filtered=dnf_filter(bates,simulated%returns,n=16,r=2)
  call assert_true(filtered%ok,'Bates filtering')

  pmd=dynamics_svm('PittMalikDoucet',phi=0.95_dp,theta=-1.0_dp,sigma=0.2_dp, &
    rho=-0.3_dp,p=0.08_dp,alpha=-0.1_dp,delta=0.15_dp)
  grids=grid_maker(pmd,n=15,k=1,r=1)
  call assert_true(all(grids%jump_counts==[0,1]),'PMD Bernoulli jump grid')
  simulated=model_simulate(pmd,80,initial_volatility=-1.0_dp,seed=123)
  filtered=dnf_filter(pmd,simulated%returns,grids=grids)
  call assert_true(filtered%ok,'PMD filtering')
  transition=transition_matrix(pmd,grids)
  call assert_close(maxval(abs(sum(transition,dim=1)-1.0_dp)),0.0_dp,1.0e-11_dp, &
    'transition columns normalize')

  dps=dynamics_svm('DuffiePanSingleton',nu=0.006_dp,omega=20.0_dp,kappa=3.0_dp, &
    theta=0.04_dp,sigma=0.3_dp,rho=-0.4_dp,rho_z=-0.5_dp)
  grids=grid_maker(dps,n=12,k=8,r=2)
  call assert_true(size(grids%jump_mid_points)==8,'DPS jump-size grid')
  call assert_true(all(grids%jump_mid_points>0.0_dp),'positive DPS jump nodes')

  write(*,'(a)') 'test_jumps: PASS'
end program test_jumps
