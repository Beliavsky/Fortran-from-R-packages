! SPDX-License-Identifier: MIT
program risk_and_pca
  use yieldcurves
  implicit none
  type(curve_t) :: curve
  type(bond_duration_result_t) :: duration
  type(zspread_result_t) :: spread
  type(pca_result_t) :: pca
  real(dp), parameter :: histories(6,3) = reshape([ &
    0.01_dp,0.02_dp,0.03_dp,0.04_dp,0.05_dp,0.06_dp, &
    0.02_dp,0.01_dp,0.04_dp,0.03_dp,0.06_dp,0.05_dp, &
    0.03_dp,0.04_dp,0.02_dp,0.05_dp,0.04_dp,0.07_dp], [6,3])

  curve = yc_curve([1.0_dp,2.0_dp,5.0_dp,10.0_dp,30.0_dp], &
    [0.040_dp,0.042_dp,0.044_dp,0.045_dp,0.046_dp])
  duration = yc_bond_duration(100.0_dp,0.05_dp,10.0_dp,0.045_dp)
  spread = yc_zspread(95.0_dp,0.04_dp,5.0_dp,curve)
  pca = yc_pca(histories,3)
  if (.not. duration%ok) error stop trim(duration%message)
  if (.not. spread%ok) error stop trim(spread%message)
  if (.not. pca%ok) error stop trim(pca%message)

  print '(a,f12.6)', 'Bond modified duration: ',duration%modified_duration
  print '(a,f12.6)', 'Z-spread:               ',spread%zspread
  print '(a,3f12.6)', 'PCA variance shares:    ',pca%variance_explained
end program risk_and_pca
