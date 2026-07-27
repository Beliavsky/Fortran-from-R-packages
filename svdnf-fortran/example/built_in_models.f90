! SPDX-License-Identifier: GPL-3.0-only
program built_in_models
  use svdnf
  implicit none
  character(len=24), parameter :: names(7) = [character(len=24) :: &
    'Heston','Bates','DuffiePanSingleton','Taylor','TaylorWithLeverage', &
    'PittMalikDoucet','CAPM_SV']
  type(svm_dynamics) :: dynamics
  type(grid_type) :: grids
  integer :: i

  do i=1,size(names)
    dynamics=dynamics_svm(trim(names(i)))
    grids=grid_maker(dynamics,n=15,k=6,r=2)
    write(*,'(a24,3i6)') trim(dynamics%model),size(grids%var_mid_points), &
      size(grids%jump_mid_points),size(grids%jump_counts)
  end do
end program built_in_models
