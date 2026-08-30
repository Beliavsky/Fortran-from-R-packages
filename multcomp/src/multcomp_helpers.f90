! SPDX-License-Identifier: GPL-2.0-only
module multcomp_helpers
  use multcomp_kinds, only : dp
  use multcomp_types, only : mtest_type, confidence_interval_type, cld_type
  use multcomp_contrasts, only : tukey_pairs
  use multcomp_cld, only : compact_letter_display
  implicit none
  private

  public :: cld_from_tukey_test
  public :: cld_from_tukey_confint

contains

  subroutine cld_from_tukey_test(test, nlevels, level, result, decreasing)
    type(mtest_type), intent(in) :: test !! Tukey pairwise test result whose p-values determine significant comparisons.
    integer, intent(in) :: nlevels !! Number of factor levels represented by the complete Tukey contrast family.
    real(dp), intent(in) :: level !! Significance cutoff; comparisons with p-value below this value are separated.
    type(cld_type), intent(out) :: result !! Compact-letter display for the factor levels.
    logical, intent(in), optional :: decreasing !! If true, reverse the final letter-column ordering convention.

    integer, allocatable :: lower_group(:)
    integer, allocatable :: upper_group(:)
    logical, allocatable :: significant(:)

    result%ok = .false.
    result%message = ''
    if (.not. test%ok) then
      result%message = 'test object is not valid'
      return
    end if
    if (level <= 0.0_dp .or. level >= 1.0_dp) then
      result%message = 'significance level must lie strictly between zero and one'
      return
    end if
    call tukey_pairs(nlevels, lower_group, upper_group)
    if (size(test%pvalue) /= size(lower_group)) then
      result%message = 'test does not contain a complete Tukey pairwise family'
      return
    end if
    significant = test%pvalue < level
    if (present(decreasing)) then
      call compact_letter_display(significant, lower_group, upper_group, nlevels, result, decreasing)
    else
      call compact_letter_display(significant, lower_group, upper_group, nlevels, result)
    end if
  end subroutine cld_from_tukey_test

  subroutine cld_from_tukey_confint(intervals, nlevels, result, decreasing)
    type(confidence_interval_type), intent(in) :: intervals !! Tukey intervals used to flag comparisons excluding zero.
    integer, intent(in) :: nlevels !! Number of factor levels represented by the complete Tukey contrast family.
    type(cld_type), intent(out) :: result !! Compact-letter display for the factor levels.
    logical, intent(in), optional :: decreasing !! If true, reverse the final letter-column ordering convention.

    integer, allocatable :: lower_group(:)
    integer, allocatable :: upper_group(:)
    logical, allocatable :: significant(:)

    result%ok = .false.
    result%message = ''
    if (.not. intervals%ok) then
      result%message = 'confidence-interval object is not valid'
      return
    end if
    call tukey_pairs(nlevels, lower_group, upper_group)
    if (size(intervals%lower) /= size(lower_group)) then
      result%message = 'confidence intervals do not contain a complete Tukey pairwise family'
      return
    end if
    significant = .not. (intervals%lower < 0.0_dp .and. intervals%upper > 0.0_dp)
    if (present(decreasing)) then
      call compact_letter_display(significant, lower_group, upper_group, nlevels, result, decreasing)
    else
      call compact_letter_display(significant, lower_group, upper_group, nlevels, result)
    end if
  end subroutine cld_from_tukey_confint

end module multcomp_helpers
