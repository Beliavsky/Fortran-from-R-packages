program test_contrasts
  use multcomp_kinds, only : dp
  use multcomp_types, only : contrast_matrix_type
  use multcomp_contrasts, only : contr_mat
  implicit none

  type(contrast_matrix_type) :: cm
  character(len=32) :: families(10)
  integer :: i
  real(dp) :: n(4)

  n = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
  families = [character(len=32) :: 'Dunnett', 'Tukey', 'Sequen', 'AVE', &
    'Changepoint', 'Williams', 'Marcus', 'McDermott', 'UmbrellaWilliams', 'GrandMean']
  do i = 1, size(families)
    call contr_mat(n, trim(families(i)), cm)
    if (.not. cm%ok) error stop 'documented contrast family construction failed'
    if (size(cm%value, 2) /= 4) error stop 'contrast family has incorrect column count'
    if (maxval(abs(sum(cm%value, dim=2))) > 2.0e-14_dp) error stop 'contrast row does not sum to zero'
  end do

  call contr_mat(n, 'Tukey', cm)
  if (size(cm%value, 1) /= 6 .or. size(cm%value, 2) /= 4) error stop 'Tukey dimensions'
  if (any(abs(cm%value(1, :) - [-1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]) > 1.0e-14_dp)) &
    error stop 'Tukey row mismatch'

  call contr_mat(n, 'AVE', cm)
  if (abs(cm%value(1, 2) + 2.0_dp / 9.0_dp) > 1.0e-12_dp) error stop 'AVE mismatch'

  call contr_mat(n, 'UmbrellaWilliams', cm)
  if (abs(cm%value(5, 2) - 0.4_dp) > 1.0e-12_dp) error stop 'UmbrellaWilliams mismatch'
  if (abs(cm%value(5, 3) - 0.6_dp) > 1.0e-12_dp) error stop 'UmbrellaWilliams mismatch'

  call contr_mat(n, 'Dunnett', cm, base=3)
  if (size(cm%value, 1) /= 3) error stop 'Dunnett base-group row count mismatch'
  if (any(abs(cm%value(:, 3) + 1.0_dp) > 1.0e-14_dp)) error stop 'Dunnett base column mismatch'
end program test_contrasts
