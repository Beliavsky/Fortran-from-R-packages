! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_compound_indicators
  use gb2_kinds, only : dp
  use gb2_compound, only : pcgb2, qcgb2, moment_cgb2, incomplete_moment_cgb2
  implicit none
  private
  public :: arpt_cgb2, arpr_cgb2, rmpg_cgb2, qsr_cgb2, main_cgb2
contains
  real(dp) function arpt_cgb2(prop,shape1,scale,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: prop,shape1,scale,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    v=prop*qcgb2(0.5_dp,shape1,scale,shape2,shape3,pl0,pl,decomp)
  end function arpt_cgb2
  real(dp) function arpr_cgb2(prop,shape1,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: prop,shape1,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    v=pcgb2(arpt_cgb2(prop,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp),shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)
  end function arpr_cgb2
  real(dp) function rmpg_cgb2(arpr,shape1,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: arpr,shape1,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    v=1.0_dp-qcgb2(arpr/2.0_dp,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)/qcgb2(arpr,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)
  end function rmpg_cgb2
  real(dp) function qsr_cgb2(shape1,shape2,shape3,pl0,pl,decomp) result(v)
    real(dp), intent(in) :: shape1,shape2,shape3,pl0(:),pl(:)
    character(len=*), intent(in), optional :: decomp
    real(dp) :: q20,q80
    q20=qcgb2(0.2_dp,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)
    q80=qcgb2(0.8_dp,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)
    v=(1.0_dp-incomplete_moment_cgb2(q80,1.0_dp,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp))/ &
      incomplete_moment_cgb2(q20,1.0_dp,shape1,1.0_dp,shape2,shape3,pl0,pl,decomp)
  end function qsr_cgb2
  subroutine main_cgb2(prop,shape1,scale,shape2,shape3,pl0,pl,values,decomp)
    real(dp), intent(in) :: prop,shape1,scale,shape2,shape3,pl0(:),pl(:)
    real(dp), intent(out) :: values(5)
    character(len=*), intent(in), optional :: decomp
    real(dp) :: ar
    ar=arpr_cgb2(prop,shape1,shape2,shape3,pl0,pl,decomp)
    values=[qcgb2(0.5_dp,shape1,scale,shape2,shape3,pl0,pl,decomp), &
      moment_cgb2(1.0_dp,shape1,scale,shape2,shape3,pl0,pl,decomp), ar, &
      rmpg_cgb2(ar,shape1,shape2,shape3,pl0,pl,decomp), qsr_cgb2(shape1,shape2,shape3,pl0,pl,decomp)]
  end subroutine main_cgb2
end module gb2_compound_indicators
