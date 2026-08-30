program test_cld_mmm
  use multcomp_kinds, only : dp
  use multcomp_types, only : cld_type, parm_type
  use multcomp_cld, only : compact_letter_display
  use multcomp_mmm, only : mmm_parm_from_iid
  use multcomp_parm, only : make_parm
  implicit none

  type(cld_type) :: letters
  type(parm_type) :: blocks(2)
  type(parm_type) :: combined
  logical :: significant(6)
  integer :: group_a(6)
  integer :: group_b(6)
  real(dp) :: iid(2, 4)
  real(dp) :: sigma1(1, 1)
  real(dp) :: sigma2(1, 1)

  group_a = [1, 1, 1, 2, 2, 3]
  group_b = [2, 3, 4, 3, 4, 4]
  significant = [.false., .true., .true., .false., .true., .false.]
  call compact_letter_display(significant, group_a, group_b, 4, letters)
  if (.not. letters%ok) error stop 'compact-letter display failed'
  call check_letter_constraints(letters%letter_matrix, significant, group_a, group_b)

  sigma1(1, 1) = 4.0_dp
  sigma2(1, 1) = 9.0_dp
  call make_parm([1.0_dp], sigma1, blocks(1))
  call make_parm([-2.0_dp], sigma2, blocks(2))
  iid(1, :) = [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp]
  iid(2, :) = [2.0_dp, -2.0_dp, 1.0_dp, -1.0_dp]
  call mmm_parm_from_iid(blocks, iid, combined)
  if (.not. combined%ok) error stop 'multiple marginal model covariance failed'
  if (maxval(abs(combined%coef - [1.0_dp, -2.0_dp])) > 1.0e-14_dp) &
    error stop 'combined coefficient mismatch'
  if (abs(combined%vcov(1, 1) - 4.0_dp) > 1.0e-13_dp) error stop 'first marginal variance changed'
  if (abs(combined%vcov(2, 2) - 9.0_dp) > 1.0e-13_dp) error stop 'second marginal variance changed'
  if (combined%vcov(1, 2) <= 0.0_dp) error stop 'expected positive cross-model covariance'

contains

  subroutine check_letter_constraints(matrix, flags, first, second)
    logical, intent(in) :: matrix(:, :) !! Group-by-letter membership matrix produced by the CLD algorithm.
    logical, intent(in) :: flags(:) !! Significance status for each tested pair.
    integer, intent(in) :: first(:) !! First group index for each tested pair.
    integer, intent(in) :: second(:) !! Second group index for each tested pair.

    integer :: i
    logical :: shares

    do i = 1, size(flags)
      shares = any(matrix(first(i), :) .and. matrix(second(i), :))
      if (flags(i) .and. shares) error stop 'significant groups share a compact-display letter'
      if (.not. flags(i) .and. .not. shares) error stop 'nonsignificant groups do not share a letter'
    end do
  end subroutine check_letter_constraints

end program test_cld_mmm
